#!/usr/bin/env bash

source environment
echo "Sourcing file \"environment\""
ACCOUNT_ID=$(aws sts get-caller-identity | jq -r .Account)

declare -A keys=( ["account_id"]=$ACCOUNT_ID ["app_name"]=$APP_NAME ["image_name"]=$IMAGE_NAME
 ["image_tag"]=$IMAGE_TAG ["dpss_task_cpu"]=$DPSS_TASK_CPU ["dpss_task_memory"]=$DPSS_TASK_MEMORY])

FILE=newjunk.tf
if [[ -f "$FILE" ]]; then
    echo "$FILE exist - overwriting it..."
fi

for val in "${!keys[@]}"
do
  printf "variable \"$val\" {\n  default = \"${keys[$val]}\" \n}\n\n" >> newjunk.tf
done

echo "Wrote Terraform configuration to file $FILE"
