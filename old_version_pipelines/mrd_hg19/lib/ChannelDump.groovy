//
// This file holds several functions used to write channel to csv file
//

class ChannelDump {

    // initialise 
    static def workflow
    static def outdir
    
    public static void initialise(workflow, outdir) {   
        this.workflow = workflow
        this.outdir = outdir
    }

    //
    // Get file full path
    //
    public static String getPath(outdir, fileName) {
        return "${outdir}/channels/${this.workflow.runName}/${fileName}"
    }

    //
    // Function to write channel to in csv file
    //

    public static void writeCsv(ChannelArray, fileHeader, fileName) {

        def outdir = this.outdir

        def filePath = getPath(outdir, fileName)

        def newFile = new File(filePath)

        if(!newFile.getParentFile().exists()) { // Checks for "full/path/to".
            newFile.getParentFile().mkdirs(); // Creates any missing directories.
        }

        newFile.text = fileHeader.join(',') + System.lineSeparator() + ChannelArray*.join(',').join(System.lineSeparator())


    }

}