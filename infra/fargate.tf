data "aws_caller_identity" "current"{}
data "aws_region" "current" {}
data "aws_vpc" "default" {
  default = true
}

data "aws_availability_zones" "available" {}

data "aws_subnet" "default" {
  count             = 3
  vpc_id            = "${data.aws_vpc.default.id}"
  availability_zone = "${data.aws_availability_zones.available.names[count.index]}"
  default_for_az    = true
}

data "aws_ecs_cluster" "default"{
  cluster_name = "default"
}

locals {
  common_tags = "${map(
    "managedBy" , "terraform",
    "Name"      , "${var.GAE_INFRA_TAG_SERVICE}-${var.GAE_INFRA_TAG_PROJECT}-${var.GAE_DEPLOYMENT_STAGE}",
    "project"   , "${var.GAE_INFRA_TAG_PROJECT}",
    "env"       , "${var.GAE_DEPLOYMENT_STAGE}",
    "service"   , "${var.GAE_INFRA_TAG_SERVICE}",
    "owner"     , "${var.GAE_INFRA_TAG_OWNER}"
  )}"
}

resource "aws_iam_role" "task-executor" {
  name = "gae-exporter-${var.GAE_DEPLOYMENT_STAGE}"
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

resource "aws_iam_role_policy_attachment" "task-executor-ecs" {
  role = "${aws_iam_role.task-executor.name}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

resource "aws_iam_role_policy" "AWS-Events-Invoke" {
  name = "gae-exporter-cloudwatch-event-invoke-${var.GAE_DEPLOYMENT_STAGE}"
  role = "${aws_iam_role.task-executor.id}"
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "ecs:RunTask"
            ],
            "Resource": [
                "arn:aws:ecs:*:${data.aws_caller_identity.current.account_id}:task-definition/gae-exporter-${var.GAE_DEPLOYMENT_STAGE}:*"
            ],
            "Condition": {
                "ArnLike": {
                    "ecs:cluster": "arn:aws:ecs:*:${data.aws_caller_identity.current.account_id}:cluster/default"
                }
            }
        },
        {
            "Effect": "Allow",
            "Action": "iam:PassRole",
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
EOF
}

resource "aws_iam_role" "task-performer" {
  name = "gae-exporter-task-performer-${var.GAE_DEPLOYMENT_STAGE}"
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

resource "aws_iam_role_policy" "task-performer" {
  name = "gae-exporter-task-performer-${var.GAE_DEPLOYMENT_STAGE}"
  role = "${aws_iam_role.task-performer.id}"
  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
            {
      "Effect": "Allow",
      "Action": "secretsmanager:Get*",
      "Resource": "arn:aws:secretsmanager:*:${data.aws_caller_identity.current.account_id}:secret:${var.GAE_SECRETS_STORE}/*"
    },
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
      "Action": [
          "dynamodb:*"
      ],
      "Resource": "arn:aws:dynamodb:*:${data.aws_caller_identity.current.account_id}:table/ga-metrics-${var.GAE_DEPLOYMENT_STAGE}"
    },
    {
      "Effect": "Allow",
      "Action":[
        "Logs:FilterLogEvents",
        "Logs:GetLogEvents",
        "Logs:GetQueryResults",
        "Logs:DescribeLogGroups",
        "Logs:DescribeLogStreams",
        "Logs:GetLogRecord",
        "Logs:StartQuery",
        "Logs:StopQuery"
      ],
      "Resource": "arn:aws:logs:*:*:*"
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

resource "aws_ecs_task_definition" "monitor" {
  # Note if changing the family, also change the IAM ROLE POLICY under ecs:RunTask; unable to reference due to dependency cycle.
  family = "gae-exporter-${var.GAE_DEPLOYMENT_STAGE}"
  execution_role_arn = "${aws_iam_role.task-executor.arn}"
  task_role_arn = "${aws_iam_role.task-performer.arn}"
  requires_compatibilities = ["FARGATE"]
  network_mode = "awsvpc"
  cpu = "256"
  memory = "512"
  container_definitions = <<DEFINITION
[
  {
    "family": "gae-exporter",
    "name": "gae-exporter-fargate",
    "image": "humancellatlas/gae-exporter-${var.GAE_DEPLOYMENT_STAGE}",
    "essential": true,
    "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-group": "${aws_cloudwatch_log_group.task-performer.name}",
          "awslogs-region": "us-east-1",
          "awslogs-stream-prefix": "ecs"
        }
    },
    "environment" :[
          {"name": "GAE_DEPLOYMENT_STAGE", "value": "${var.GAE_DEPLOYMENT_STAGE}"},
          {"name": "GAE_SECRETS_STORE", "value": "${var.GAE_SECRETS_STORE}"}
    ]
  }
]
DEFINITION
  tags = "${local.common_tags}"
}


resource "aws_cloudwatch_log_group" "task-performer" {
  name              = "/aws/service/gae-exporter-task-performer-${var.GAE_DEPLOYMENT_STAGE}"
  retention_in_days = 1827
}

resource "aws_cloudwatch_event_rule" "gae-exporter" {
  name = "gae-exporter-trigger-${var.GAE_DEPLOYMENT_STAGE}"
  description = "schedule trigger for data input to GAE"
  tags = "${local.common_tags}"
  schedule_expression = "cron(1/10 * * * ? *)"
}


resource "aws_cloudwatch_event_target" "scheduled_task" {
  rule       = "${aws_cloudwatch_event_rule.gae-exporter.name}"
  arn        = "${data.aws_ecs_cluster.default.arn}"
  role_arn   = "${aws_iam_role.task-executor.arn}"
  ecs_target = {
    task_count          = 1
    task_definition_arn = "${aws_ecs_task_definition.monitor.arn}"
    launch_type         = "FARGATE"
    platform_version    = "LATEST"

    network_configuration = {
      assign_public_ip = true
      subnets          = ["${data.aws_subnet.default.*.id}"]
    }
  }
}