terraform {
  required_version = ">=0.12"
}

variable "app_name" {
  default = "data-portal-summary-stats"
}

variable "image_name" {
  default = "data-portal-summary-stats"
}

provider "aws" {
  assume_role {
    role_arn = var.role_arn
  }
  region = var.region
}

data "aws_caller_identity" "current"{}
data "aws_region" "current" {}

// Fetch AZs in current region.
data "aws_availability_zones" "available" {}

// Virtual Private Cloud and subnets
data "aws_vpc" "data-portal-summary-stats" {
  id = var.dpss_vpc_id
}

data "aws_subnet" "dpss-sn" {
  vpc_id = "${data.aws_vpc.data-portal-summary-stats.id}"
  cidr_block = "10.0.1.0/24"
}

resource "aws_ecs_cluster" "dpss-cluster"{
  name = var.cluster_name
}

locals {
  common_tags = "${map(
    "managedBy" , "terraform"
  )}"
}

/*
In the following we first define an IAM role, then we create policies, and finally, we attach
those the policies to the roles.
*/

/*
Role and policy for ECS events.
*/
resource "aws_iam_role" "data-portal-summary-stats-ecs-events" {
  name = "data-portal-summary-stats_ecs_events"
  description = "Run dpss ecs task"
  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": [
          "ecs.amazonaws.com",
          "ecs-tasks.amazonaws.com",
          "events.amazonaws.com"
        ]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF
}

resource "aws_iam_policy" "data-portal-summary-stats-ecs-events-policy" {
  name = "data-portal-summary-stats-policy-run-task"
  description = "Run dpss ecs task"

  policy = <<DOC
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "ecs:RunTask"
            ],
            "Resource": [
                "arn:aws:ecs::${data.aws_caller_identity.current.account_id}:task-definition/data-portal-summary-stats-${var.deployment_stage}:*"
            ],
            "Condition": {
                "ArnLike": {
                     "ecs:cluster": "${aws_ecs_cluster.dpss-cluster.arn}"
                }
            }
        },
        {
            "Effect": "Allow",
            "Action": [
                "iam:PassRole"
            ],
            "Resource": [
                "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/ecsTaskExecutionRole"
            ],
            "Condition": {
                "StringLike": {
                    "iam:PassedToService": "ecs-tasks.amazonaws.com"
                }
            }
        }
    ]
}
DOC
}

// Connect role to policy.
resource "aws_iam_role_policy_attachment" "ecs-events-attach" {
  policy_arn = "${aws_iam_policy.data-portal-summary-stats-ecs-events-policy.arn}"
  role       = "${aws_iam_role.data-portal-summary-stats-ecs-events.id}"
}

resource "aws_iam_role_policy_attachment" "ecs-events-attach2" {
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
  role       = "${aws_iam_role.data-portal-summary-stats-ecs-events.id}"
}

// This attachment is required to pull the container from ECR.
resource "aws_iam_role_policy_attachment" "ecs-events-attach3" {
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceEventsRole"
  role       = "${aws_iam_role.data-portal-summary-stats-ecs-events.id}"
}

/*
Role and policies for the task performer.
*/
resource "aws_iam_role" "data-portal-summary-stats-task-performer" {
  name = "data-portal-summary-stats-task-performer"
  tags = "${local.common_tags}"
  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": [
          "ecs-tasks.amazonaws.com"
        ]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF
}

resource "aws_iam_policy" "data-portal-summary-stats-task-performer-policy" {
  name        = "data-portal-summary-stats-task-performer-policy"
  description = "Perform task"

  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "tag:GetTagKeys",
        "tag:GetResources",
        "tag:GetTagValues",
        "cloudwatch:*"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action":[
        "logs:FilterLogEvents",
        "logs:GetLogEvents",
        "logs:GetQueryResults",
        "logs:DescribeLogGroups",
        "logs:DescribeLogStreams",
        "logs:GetLogRecord",
        "logs:StartQuery",
        "logs:StopQuery"
      ],
      "Resource": "arn:aws:logs:*:*:*"
    },
    {
      "Effect": "Allow" ,
      "Action": [
        "ecs:CreateCluster",
        "ecs:DeregisterContainerInstance",
        "ecs:DiscoverPollEndpoint",
        "ecs:Poll",
        "ecs:RegisterContainerInstance",
        "ecs:StartTelemetrySession",
        "ecs:Submit*",
        "ecs:ListTasks"
      ],
      "Resource": "*"
    },
    {
      "Sid": "VisualEditor0",
      "Effect": "Allow",
      "Action": "s3:ListBucket",
      "Resource": "arn:aws:s3:::*"
    },
    {
      "Sid": "VisualEditor1",
      "Effect": "Allow",
      "Action": [
          "s3:PutObject",
          "s3:GetObject"
      ],
      "Resource": "arn:aws:s3:::*/*"
    }
  ]
}
EOF
}

// Attached policies to role.
resource "aws_iam_role_policy_attachment" "task-performer-attach" {
  policy_arn = "${aws_iam_policy.data-portal-summary-stats-task-performer-policy.arn}"
  role       = "${aws_iam_role.data-portal-summary-stats-task-performer.id}"
}

resource "aws_ecs_task_definition" "dpss_ecs_task_definition" {
  family                   = "${var.app_name}-${var.deployment_stage}"
  execution_role_arn       = "${aws_iam_role.data-portal-summary-stats-ecs-events.arn}"
  task_role_arn            = "${aws_iam_role.data-portal-summary-stats-task-performer.arn}"
  requires_compatibilities = ["FARGATE"]
  network_mode             = "awsvpc"
  cpu                      = var.dpss_task_cpu
  memory                   = var.dpss_task_memory
  tags                     = "${local.common_tags}"

  container_definitions    = <<DEFINITION
[
  {
    "family": "${var.app_name}",
    "name": "${var.app_name}",
    "image": "${var.ecr_path}${var.image_name}:${var.image_tag}",
    "essential": true,
    "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-group": "${aws_cloudwatch_log_group.task-execution.name}",
          "awslogs-region": "${var.region}",
          "awslogs-stream-prefix": "ecs"
        }
    },
    "command": [
          "--environ",
          "dev",
          "--source",
          "canned",
          "--blacklist",
          "false",
          "--min_gene_count",
          "1800"
     ]
  }
]
DEFINITION
}

// To run ECS scheduled tasks we need to use CloudWatch event rules...
// "cron(1/2 * * * ? *)": every 2 min
// "cron(0 0 * * ? *)": every day at midnight (region's timezone)
// "rate(6 hours)": every 6 h, starting from invocation
resource "aws_cloudwatch_event_rule" "dpss-scheduler" {
  name                = "dpss-trigger-${var.deployment_stage}"
  description         = "Schedule to run data-portal-summary-stats"
  tags                = "${local.common_tags}"
  schedule_expression = "cron(30 2 * * ? *)"
}

resource "aws_cloudwatch_event_target" "scheduled_task" {
  target_id = "run-scheduled-dpss-task-every-24h"
  arn       = "${aws_ecs_cluster.dpss-cluster.arn}"
  rule      = "${aws_cloudwatch_event_rule.dpss-scheduler.name}"
  role_arn  = "${aws_iam_role.data-portal-summary-stats-ecs-events.arn}"

  ecs_target {
      task_count          = 1
      task_definition_arn = "${aws_ecs_task_definition.dpss_ecs_task_definition.arn}"
      launch_type         = "FARGATE"
      platform_version    = "LATEST"

      network_configuration {
        assign_public_ip  = true
        subnets           = "${data.aws_subnet.dpss-sn.*.id}"
        security_groups   = ["${var.dpss_security_group_id}"]
      }
  }
  input = <<DOC
{
  "containerOverrides": [
    {
      "name": "${var.app_name}",
      "command": [
        "--environ", "${var.deployment_stage}",
        "--source","${var.matrix_source}",
        "--blacklist","${var.blacklist}",
        "--min_gene_count","${var.min_gene_count}"
      ]
    }
  ]
}
DOC
}

resource "aws_cloudwatch_log_group" "task-execution" {
  name              = "/ecs/${var.app_name}-${var.deployment_stage}"
  retention_in_days = 1827  // that's 5 years
}
