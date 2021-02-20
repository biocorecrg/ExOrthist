aws ssm send-command \
	--document-name "AWS-RunShellScript" \
	--targets "Key=InstanceIds,Values=instance-id" \
	--cli-input-json file://run.json
