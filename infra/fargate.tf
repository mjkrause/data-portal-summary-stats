terraform {
  required_version = ">=0.12"
}

provider "aws" {
  assume_role {
    role_arn = var.role_arn
  }
  region = var.region
}

data "aws_caller_identity" "current"{}
data "aws_region" "current" {}

//data "aws_vpc" "default" {
//  default = true
//}

data "aws_vpc" "data-portal-summary-stats" {
  id = var.dpss_vpc_id
}

data "aws_iam_role" "ecs_task_execution_role" {
  name = "ecsTaskExecutionRole"
}

// Fetch AZs in current region.
data "aws_availability_zones" "available" {}

data "aws_subnet" "default" {
  vpc_id = "${data.aws_vpc.data-portal-summary-stats.id}"
  cidr_block = "10.0.1.0/24"
}

//data "aws_subnet" "default" {
//  count = 2
//  vpc_id = "${data.aws_vpc.default.id}"
//  availability_zone = "${data.aws_availability_zones.available.names[count.index]}"
//  default_for_az = true
//}

data "aws_ecs_cluster" "default"{
  cluster_name = var.app_name
}

locals {
  common_tags = "${map(
    "managedBy"       , "terraform",
    "environ"         , "dev",
    "source"          , "canned",
    "min_gene_count"  , "1200",
    "blacklist"       , "true"
  )}"
}

/*
In the following we first define an IAM role, then we create two policies, one for ECR/ECS
service, another for S3 service. Finally, we attached those two policies to the role.
*/

resource "aws_iam_role" "data-portal-summary-stats-role-tf" {
  name = "data-portal-summary-stats-role-tf"
  description = "Allows ECS tasks to call AWS services on your behalf."
  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": "sts:AssumeRole",
      "Effect": "Allow",
      "Principal": {
        "Service": [
          "ecs.amazonaws.com",
          "ecs-tasks.amazonaws.com",
          "s3.amazonaws.com"
        ]
      },
      "Sid": ""
    }
  ]
}
EOF
}

resource "aws_iam_policy" "data-portal-summary-stats-policy-ecs-tf" {
  name = "data-portal-summary-stats-policy-ecs-tf"
  description = "data-portal-summary-stats ECS policy (TF)"
  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "ecs:CreateCluster",
        "ecs:DeregisterContainerInstance",
        "ecs:DiscoverPollEndpoint",
        "ecs:Poll",
        "ecs:RegisterContainerInstance",
        "ecs:StartTelemetrySession",
        "ecs:Submit*",
        "ecr:GetAuthorizationToken",
        "ecr:BatchCheckLayerAvailability",
        "ecr:GetDownloadUrlForLayer",
        "ecr:BatchGetImage",
        "logs:CreateLogStream",
        "logs:PutLogEvents"
      ],
      "Resource": "*"
    }
  ]
}
EOF
}

resource "aws_iam_policy" "data-portal-summary-stats-policy-s3-tf" {
  name        = "data-portal-summary-stats-policy-s3-tf"
  description = "data-portal-summary-stats S3 policy (TF)"
  policy      = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
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


resource "aws_iam_role_policy_attachment" "dpss-attach1" {
  policy_arn = "${aws_iam_policy.data-portal-summary-stats-policy-ecs-tf.arn}"
  role       = "${aws_iam_role.data-portal-summary-stats-role-tf.name}"
}

resource "aws_iam_role_policy_attachment" "dpss-attach2" {
  policy_arn = "${aws_iam_policy.data-portal-summary-stats-policy-s3-tf.arn}"
  role       = "${aws_iam_role.data-portal-summary-stats-role-tf.name}"
}


/*
Policy and role for ECS events
*/
resource "aws_iam_role" "data-portal-summary-stats_ecs_events" {
  name = "data-portal-summary-stats_ecs_events"
  description = "Run dpss ecs task"
  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "",
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

resource "aws_iam_role_policy" "dpss_ecs-policy-run-task" {
  name = "data-portal-summary-stats-policy-run-task"
  //description = "Run dpss ecs task"
  role = "${aws_iam_role.data-portal-summary-stats_ecs_events.id}"
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
                  "${replace(aws_ecs_task_definition.dpss_ecs_task_definition.arn, "/:\\d+$/", ":*")}"
            ],
            "Condition": {
                "ArnLike": {
                    "ecs:cluster": "arn:aws:ecs:*:${var.role_arn}:cluster/data-portal-summary-stats"
                }
            }
        },
        {
            "Effect": "Allow",
            "Action": [
                "iam:ListInstanceProfiles",
                "iam:ListRoles",
                "iam:PassRole"
            ],
            "Resource": [
                "*"
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

//// Connect role to policy.
//resource "aws_iam_role_policy_attachment" "ecs_events_attach" {
//  policy_arn = "${aws_iam_policy.dpss_ecs-policy-run-task.arn}"
//  //role = "${aws_iam_role.ecs_events.name}"
//  role = "${aws_iam_role.data-portal-summary-stats_ecs_events.id}"
//}

/*
??
*/
resource "aws_iam_role" "data-portal-summary-stats-task-performer" {
  name = "${var.app_name}-${var.deployment_stage}"
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

resource "aws_iam_role_policy" "data-portal-summary-stats-task-performer" {
  name = "${var.app_name}-${var.deployment_stage}"
  role = "${aws_iam_role.data-portal-summary-stats-task-performer.id}"
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
      "Effect": "Allow" ,
      "Action": [
        "ecs:ListTasks"
      ],
      "Resource": "*"
    }
  ]
}
EOF
}

resource "aws_ecs_task_definition" "dpss_ecs_task_definition" {
  family = "data-portal-summary-stats-${var.deployment_stage}"
  execution_role_arn       = "${data.aws_iam_role.ecs_task_execution_role.arn}"
  task_role_arn            = "${aws_iam_role.data-portal-summary-stats-role-tf.arn}"
  requires_compatibilities = ["FARGATE"]
  network_mode             = "awsvpc"
  cpu                      = "2048"
  memory                   = "16384"
  tags                     = "${local.common_tags}"
  container_definitions    = <<DEFINITION
[
  {
    "family": "data-portal-summary-stats",
    "name": "data-portal-summary-stats-fargate",
    "image": "${var.ecr_path}${var.image_name}:${var.image_tag}",
    "essential": true,
    "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-group": "/ecs/${var.app_name}-${var.deployment_stage}",
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
          "1200"
     ],
    "name": "${var.app_name}"
  }
]
DEFINITION
}

// To run ECS scheduled tasks we need to use CloudWatch event rules...
resource "aws_cloudwatch_event_rule" "dpss-scheduler" {
  name                = "dpss-trigger-${var.deployment_stage}"
  description         = "Schedule to run data-portal-summary-stats"
  tags                = "${local.common_tags}"
  schedule_expression = "cron(1/5 * * * ? *)"
}

resource "aws_cloudwatch_event_target" "scheduled_task" {
  target_id = "run-scheduled-dpss-task-every-24h"
  arn = "${data.aws_ecs_cluster.default.arn}"
  rule = "${aws_cloudwatch_event_rule.dpss-scheduler.name}"
  role_arn = "${aws_iam_role.data-portal-summary-stats_ecs_events.arn}"

  ecs_target {
      task_count          = 1
      task_definition_arn = "${aws_ecs_task_definition.dpss_ecs_task_definition.arn}"
      launch_type         = "FARGATE"
      platform_version    = "LATEST"

      network_configuration {
        assign_public_ip  = true
        subnets           = "${data.aws_subnet.default.*.id}"
        security_groups = ["${var.dpss_security_group_id}"]
      }
    }
  input = <<DOC
{
  "containerOverrides": [
    {
      "command": ["--environ","dev",
                  "--source","canned",
                  "--blacklist","false",
                  "--min_gene_count","300"]
    }
  ]
}
DOC
}

resource "aws_cloudwatch_log_group" "task-execution" {
  name              = "/ecs/${var.app_name}-${var.deployment_stage}"
  retention_in_days = 1827
}
