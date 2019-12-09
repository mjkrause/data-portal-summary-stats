#!/bin/bash

export DPSS_IMAGE_TAG=0.9.3

export DPSS_DEPLOYMENT_STAGE=dev

export DPSS_MATRIX_SOURCE=fresh

export DPSS_BLACKLIST=1

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
export TF_VAR_image_tag=$DPSS_IMAGE_TAG
export TF_VAR_deployment_stage=$DSS_DEPLOYMENT_STAGE
export TF_VAR_matrix_source=$DPSS_MATRIX_SOURCE
export TF_VAR_blacklist=$DPSS_BLACKLIST
export TF_VAR_cluster_name=$DPSS_CLUSTER_NAME
export TF_VAR_dpss_task_cpu=$DPSS_TASK_CPU
export TF_VAR_dpss_task_memory=$DPSS_TASK_MEMORY
export TF_VAR_dpss_security_group_id=$DPSS_SECURITY_GROUP_ID
export TF_VAR_dpss_vpc_id=$DPSS_VPC_ID
