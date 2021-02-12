Example recipe for running an AWS Batch with CPU node instances

* Entrypoint AMI can be created following instructions here: https://www.nextflow.io/docs/latest/awscloud.html#custom-ami

export TF_VAR_amientrypoint=ami-xxxxxx
export TF_VAR_amibatch=ami-xxxxx

== Instructions ==

Go to the directory with the customized Terraform scripts and check they are valid:

```terraform validate```

You may need to ```terraform init``` first.



== Additional links ==
* [https://blog.gruntwork.io/a-comprehensive-guide-to-managing-secrets-in-your-terraform-code-1d586955ace1 Storing sensitive data with Terraform]
