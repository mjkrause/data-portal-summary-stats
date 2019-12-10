#!/bin/bash


set -o allexport
source environment.env
set +o allexport

export DPSS_CLUSTER_NAME=data-portal-summary-stats-fargate

export DPSS_TASK_CPU=2048

export DPSS_TASK_MEMORY=16384

export DPSS_SECURITY_GROUP_ID=sg-095dc7a781d1d6744

export DPSS_VPC_ID=vpc-0b94af0287c8aff49

# Make src readable by tests
export PYTHONPATH="./src:./test"

# Set environment variables for Terraform:
export TF_VAR_acc_number=$(aws sts get-caller-identity | jq -r .Account)
export TF_VAR_aws_region=$(aws configure get region)
export TF_VAR_dpss_deployment_stage=$DPSS_MTX_TARGET_STAGE
export TF_VAR_image_tag=$DPSS_IMAGE_TAG
export TF_VAR_cluster_name=$DPSS_CLUSTER_NAME
export TF_VAR_dpss_task_cpu=$DPSS_TASK_CPU
export TF_VAR_dpss_task_memory=$DPSS_TASK_MEMORY
export TF_VAR_dpss_security_group_id=$DPSS_SECURITY_GROUP_ID
export TF_VAR_dpss_vpc_id=$DPSS_VPC_ID
