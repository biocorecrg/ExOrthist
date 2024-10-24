// Validation, initialization and completion of the pipeline
// Bits from:
// nf-core https://github.com/nf-core/demo
// BioNextflow https://github.com/biocorecrg/BioNextflow

include { paramsHelp; paramsSummaryLog; paramsSummaryMap; validateParameters } from 'plugin/nf-schema'


workflow PIPELINE_INITIALISATION {

    take:
    params
    args

    main:

    log_main = """
╔╦╗┬ ┬┌─┐  ╔═╗─┐ ┬╔═╗┬─┐┌┬┐┬ ┬┬┌─┐┌┬┐
 ║ ├─┤├┤   ║╣ ┌┴┬┘║ ║├┬┘ │ ├─┤│└─┐ │
 ╩ ┴ ┴└─┘  ╚═╝┴ └─╚═╝┴└─ ┴ ┴ ┴┴└─┘ ┴

==============================================================================
annotations (GTF files)             : ${params.annotations}
genomes (fasta files)               : ${params.genomes}
cluster file (txt files)            : ${params.cluster}
pairwise evo distances              : ${params.evodists}
long distance parameters            : ${params.long_dist}
medium distance parameters          : ${params.medium_dist}
short distance parameters           : ${params.short_dist}
pre-computed alignments             : ${params.prevaln}
alignment number                    : ${params.alignmentnum}
orthogroup number                   : ${params.orthogroupnum}
extraexons (e.g. from VastDB)       : ${params.extraexons}
bona fide orthologous exon pairs    : ${params.bonafide_pairs}
orthopairs                          : ${params.orthopairs}
output (output folder)              : ${params.output}
email for notification              : ${params.email}
hook_url                            : ${params.hook_url}

INFORMATION ABOUT OPTIONS:
The long, medium, short distance cut-offs are in the format: "int_num;ex_seq;ex_len;prot_sim".
Only exon matches respecting all cut-offs are considered homologous.
- int_num (0,1,2): Number of surrounding intron positions required to be conserved.
- ex_seq (from 0 to 1): Minimum sequence similarity % between a
     pair of homologous exons and their corresponding upstream and
     downstream exons.
- ex_len (from 0 to 1): Maximum size difference between two homologous exons
     (as a fraction of either exon).
- prot_sim (from 0 to 1): Minimum sequence similarity over the entire pairwise alignment
     for a pair of protein isoforms to be considered for comparison.

See online README at https://github.com/biocorecrg/ExOrthist for further information about the options.
"""

log_plot = """
Executing with the following parameters:

output main                         : ${params.output}
output plot                         : ${params.output_plot}
geneID                              : ${params.geneID}
isoformID                           : ${params.isoformID}
relevant exons                      : ${params.relevant_exs}
reclustered gene orthology file     : ${params.sub_orthologs}
email for notification              : ${params.email}
hook_url                            : ${params.hook_url}
"""

    if (params.validate_params) {
        validateParameters()
    }

    // TODO: Consider log.info paramsSummaryLog(workflow)

    if (params.wf == "plot") {
        println log_plot
    } else {
        println log_main
    }
}

workflow PIPELINE_COMPLETION {

    take:
    subworkflow
    email
    hook_url

    main:

    workflow.onComplete {
        def text = final_message("ExOrthist", subworkflow)
        println text
        if (email) {
            sendMail(to: email, subject: "[ExOrthist] Execution finished", body: text)
        }
        if (hook_url) {
            notify_slack(text, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please check documentation or file an issue at https://github.com/biocorecrg/ExOrthist/issues"
    }
}

def trim_NF_date(nfdate){
    newdate = nfdate.substring(0, nfdate.indexOf(".")).replaceAll("T", " ")
    return(newdate)
}

// notify
def notify_slack(text, hook) {
    def myFile = file('./notify.json')
    myFile << '{"text": "'
    myFile << text.replace("\n",'\\n')
    myFile << '"}'
    println "curl -X POST -H 'Content-type: application/json' -d @./notify.json ${hook}".execute().text
    myFile.delete()

}

def final_message(title="", subworkflow="") {
	def ostart = "${workflow.start}"
	def ostop = "${workflow.complete}"
	def start = trim_NF_date(ostart)
        def stop = trim_NF_date(ostop)
        def error = ""
        if (workflow.errorReport) {
            error = "\n```${workflow.errorReport}```\n"
        }
        def extrasw = ""
        if (subworkflow) {
            extrasw = " ($subworkflow)"
        }

	def message =  "-"*51 + "\n"
	message = message + "*Pipeline ${title}${extrasw} completed!*".center(51) + "\n"
        message = message + "-"*51 + "\n"

	message = message + "- Launched by `$workflow.userName`" + "\n"
	message = message + "- Started at $start" + "\n"
    	message = message + "- Finished at $stop" + "\n"
    	message = message + "- Time elapsed: $workflow.duration" + "\n"
    	message = message + "- Execution status: ${ workflow.success ? 'OK' : 'failed' }" + "\n"
    	message = message + "```$workflow.commandLine```"+ "\n"
    	message = message + error + "-"*51 + "\n"
return (message)

}
