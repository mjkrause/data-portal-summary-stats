#!/bin/bash

source environment.env

# Make src readable by tests
export PYTHONPATH="./src:./test"

# Set environment variables for Terraform:
export TF_VAR_acc_number=$(aws sts get-caller-identity | jq -r .Account)
export TF_VAR_aws_region=$(aws configure get region)
export TF_VAR_image_tag=$IMAGE_TAG
export TF_VAR_deployment_stage=$DSS_DEPLOYMENT_STAGE
export TF_VAR_matrix_source=$DPSS_MATRIX_SOURCE
export TF_VAR_blacklist=$DPSS_BLACKLIST
export TF_VAR_cluster_name=$DPSS_CLUSTER_NAME
export TF_VAR_dpss_task_cpu=$DPSS_TASK_CPU
export TF_VAR_dpss_task_memory=$DPSS_TASK_MEMORY
export TF_VAR_dpss_security_group_id=$DPSS_SECURITY_GROUP_ID
export TF_VAR_dpss_vpc_id=$DPSS_VPC_ID
