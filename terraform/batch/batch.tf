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


resource "aws_batch_compute_environment" "nf-spot" {

    compute_environment_name = "nf-spot"

    compute_resources {
      // instance_role = "${aws_iam_instance_profile.S3access.arn}"
      instance_role = "arn:aws:iam::132458770246:instance-profile/S3access"
      bid_percentage = var.bid_percentage

      image_id = var.amibatch

      max_vcpus = var.max_vcpus
      min_vcpus = var.min_vcpus
      desired_vcpus = var.desired_vcpus

      instance_type = [ "optimal" ]

      subnets = ["subnet-8a280df7", "subnet-c54d6588", "subnet-b85ab5d2"]

      spot_iam_fleet_role = "arn:aws:iam::132458770246:role/AmazonEC2SpotFleetRole"

      type = "SPOT"

      security_group_ids = [ aws_security_group.allow_all.id ]

    }

    // service_role = "${aws_iam_role.AWSBatchServiceRole.arn}"
    service_role = "arn:aws:iam::132458770246:role/service-role/AWSBatchServiceRole"
    type         = "MANAGED"
    depends_on   = [ aws_iam_policy_attachment.AWSBatchServiceRole-policy-attachment ]


}

resource "aws_batch_job_queue" "spot" {
  name                 = "spot"
  state                = "ENABLED"
  priority             = 1
  compute_environments = [ aws_batch_compute_environment.nf-spot.arn ]

}
