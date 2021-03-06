//batch.tf

variable "amibatch" {
  type = string
}

variable "bid_percentage" {
  type = number
  default = 50
}

variable "max_vcpus" {
  type = number
  default = 16
}

variable "min_vcpus" {
  type = number
  default = 0
}

variable "desired_vcpus" {
  type = number
  default = 0
}

variable "instance_batch" {
  type    = list(string)
  default = ["optimal"]
}

resource "aws_batch_compute_environment" "nf-compute-spot" {

    compute_environment_name = "nf-compute-spot"

    compute_resources {
      instance_role = aws_iam_instance_profile.ComputeInstanceProfile.arn
      bid_percentage = var.bid_percentage

      image_id = var.amibatch

      max_vcpus = var.max_vcpus
      min_vcpus = var.min_vcpus
      desired_vcpus = var.desired_vcpus

      instance_type = var.instance_batch

      // TODO: Eventually migrate subnets as well
      subnets = ["subnet-8a280df7", "subnet-c54d6588", "subnet-b85ab5d2"]

      spot_iam_fleet_role = aws_iam_role.ClusterFleetRole.arn

      type = "SPOT"

      security_group_ids = [ aws_security_group.allow_all.id ]

    }

    service_role = aws_iam_role.ClusterRole.arn
    type         = "MANAGED"
    depends_on   = [aws_iam_policy_attachment.AWSBatchServiceRole-policy-attachment]

    tags = {
  	 name = "nf-cluster"
    }
}

resource "aws_batch_job_queue" "spot" {

  name                 = "spot"
  state                = "ENABLED"
  priority             = 1
  compute_environments = [ aws_batch_compute_environment.nf-compute-spot.arn ]

  tags = {
   name = "nf-queue"
  }
}
