// From https://github.com/biocorecrg/BioNextflow/

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

