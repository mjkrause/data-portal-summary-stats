#!/usr/bin/env bash

set -euo pipefail  # -e: exit on first failed command

# Check if AWS profile has been set:
VAR1="$(head -n 1 environment)"
VAR2="AWS_PROFILE=\"\""  # AWS profile not set
if [ "$VAR1" = "$VAR2" ]; then
  echo "Aborting: AWS_PROFILE in file environment cannot be empty"
  exit 1
fi

source environment

ACCOUNT_ID=$(aws sts get-caller-identity | jq -r .Account)
AWS_REGION=$(aws configure get region)

declare -A keys=(
  ["acc_number"]=$ACCOUNT_ID
  ["region"]=$AWS_REGION
  ["app_name"]=$APP_NAME
  ["role_arn"]="arn:aws:iam::${ACCOUNT_ID}:role/dcp-developer"
  ["ecr_path"]="${ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/"
  ["image_name"]=$IMAGE_NAME
  ["image_tag"]=$IMAGE_TAG
  ["deployment_stage"]=$DEPLOYMENT_STAGE
  ["cluster_name"]=$CLUSTER_NAME
  ["dpss_task_cpu"]=$DPSS_TASK_CPU
  ["dpss_task_memory"]=$DPSS_TASK_MEMORY
)

VARS_TF=./infra/variables.tf
if [[ -f "$VARS_TF" ]]; then
    echo "$VARS_TF exist - overwriting it..."
    rm $VARS_TF
fi

for val in "${!keys[@]}"
do
  # ToDo: check if val is empty
  printf "variable \"$val\" {\n  default = \"${keys[$val]}\" \n}\n\n" >> $VARS_TF
done

echo "Wrote Terraform configuration to file $VARS_TF"
