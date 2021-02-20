Example recipe for running an AWS Batch with CPU node instances

* Entrypoint AMI can be created following instructions here: https://www.nextflow.io/docs/latest/awscloud.html#custom-ami

Specific variables are stored in BASH variables, such as:

export TF_VAR_amientrypoint=ami-xxxxxx
export TF_VAR_amibatch=ami-xxxxx

A good practice is keeping it in a file (e.g. <code>myvars.env</code>) and declare them by ```source myvars.env``` before running the commands below.

== Instructions ==

Go to the directory with the customized Terraform scripts and check they are valid:

```terraform validate```

You may need to ```terraform init``` first.

To review which infrastructure will be deployed:

```terraform plan```

To actually deploy it:

```terraform apply```

--- BELOW TO ADD System manager AWS CLI step --- 


Once work is done:

```terraform destroy```


== Additional links ==
* [https://blog.gruntwork.io/a-comprehensive-guide-to-managing-secrets-in-your-terraform-code-1d586955ace1 Storing sensitive data with Terraform]
