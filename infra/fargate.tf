data "aws_caller_identity" "current"{}
data "aws_region" "current" {}
data "aws_vpc" "default" {
  default = true
}

data "aws_availability_zones" "available" {}

data "aws_subnet" "default" {
  count             = 2
  vpc_id            = "${data.aws_vpc.default.id}"
  availability_zone = "${data.aws_availability_zones.available.names[count.index]}"
  default_for_az    = true
}

data "aws_ecs_cluster" "default"{
  cluster_name = "default"
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

resource "aws_iam_role" "data-portal-summary-stats-role" {
  name = "data-portal-summary-stats-role"
  assume_role_policy = <<EOF
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
                "ecr:BatchGetImage"
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

resource "aws_ecs_task_definition" "monitor" {
  family = "data-portal-summary-stats-${var.deployment_stage}"
  execution_role_arn = var.execution_role_arn
  task_role_arn = var.task_role_arn
  requires_compatibilities = ["FARGATE"]
  network_mode = "awsvpc"
  cpu = "2048"
  memory = "16384"
  revision = 26
  container_definitions = <<DEFINITION
[
  {
    "family": "data-portal-summary-stats",
    "name": "data-portal-summary-stats-fargate",
    "image": "${var.arn_}${var.image_name}:${var.image_tag}"
    "essential": true,
    "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-group": "/ecs/${var.resource_name}",
          "awslogs-region": "${var.region}"",
          "awslogs-stream-prefix": "ecs"
        }
    },
    "command": [
          "--environ",
          "prod",
          "--source",
          "fresh",
          "--blacklist",
          "true",
          "--min_gene_count",
          "1200"
     ],
    "name": "${var.resource_name}"
  }
]
DEFINITION
  tags = "${local.common_tags}"
}

resource "aws_cloudwatch_log_group" "task-execution" {
  name              = "/ecs/${var.resource_name}-${var.deployment_stage}"
  retention_in_days = 1827
}

resource "aws_cloudwatch_event_rule" "dpss-scheduler" {
  name = "dpss-trigger-${var.deployment_stage}"
  description = "Schedule to run data-portal-summary-stats"
  tags = "${local.common_tags}"
  schedule_expression = "cron(*/5 * * * *)"
}

resource "aws_cloudwatch_event_target" "scheduled_task" {
  rule       = "${aws_cloudwatch_event_rule.dpss-scheduler.name}"
  arn        = "${data.aws_ecs_cluster.default.arn}"
  role_arn   = "${aws_iam_role.data-portal-summary-stats-role.arn}"

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

resource "aws_cloudwatch_event_target" "" {
  arn = ""
  rule = "${aws_cloudwatch_event_rule.dpss-scheduler.name}"
}