//
// This file holds several functions used to perform input csv validation.
//

import org.everit.json.schema.Schema
import org.everit.json.schema.loader.SchemaLoader
import org.everit.json.schema.ValidationException
import org.json.JSONObject
import org.json.JSONTokener
import org.json.JSONArray
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

class CsvSchema {

    // initialise workflow and log 
    static def workflow
    static def log
    
    public static void initialise(workflow, log) {   
        this.workflow = workflow
        this.log = log
    }

    //
    // Resolve Schema path relative to main workflow directory
    //
    public static String getSchemaPath(workflow, schema_filename) {
        return "${workflow.projectDir}/${schema_filename}"
    }

    //
    // Function to check all columns in csv file
    //

    public static void validateCsv(csvArray, schema_filename) {

        def workflow = this.workflow
        def log = this.log

        def has_error = false
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Collect expected parameters from the schema

        def json = new File(getSchemaPath(workflow, schema_filename=schema_filename)).text

        def Map schemaColumns = (Map) new JsonSlurper().parseText(json).get('items')

        def enums = [:]

        for (p in schemaColumns.get('properties')) {
            if (schemaColumns.get('properties')[p.key].containsKey('enum')) {
                enums[p.key] = schemaColumns.get('properties')[p.key]['enum']
            }
        }


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Convert csvArray
        def newCsvArray = []
        for (row in csvArray) {
            def tmpMap = [:]
            row.each { key, cell ->
                def value
                if (cell == null) {
                    value = cell
                } else if (cell ==~ /^\d+$/){
                    value = cell.toInteger()
                } else if (cell ==~ /[+-]?[0-9]*[.]?[0-9]+/) {
                    value = cell.toFloat()
                } else if (cell.toLowerCase() == "true" ) {
                    value = true
                } else if (cell.toLowerCase() == "false" ) {
                    value = false
                } else {
                    value = cell
                }
                if (value != null) {
                    tmpMap[key] = value
                }
            }
            newCsvArray.add(tmpMap)
        }
            
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Validate parameters against the schema
        InputStream input_stream = new File(getSchemaPath(workflow, schema_filename=schema_filename)).newInputStream()
        JSONObject raw_schema = new JSONObject(new JSONTokener(input_stream))

        Schema schema = SchemaLoader.load(raw_schema)

        // Convert to JSONObject
        def jsonCsv = new JsonBuilder(newCsvArray)
        JSONArray csv_json = new JSONArray(jsonCsv.toString())

        // Validate
        try {
            schema.validate(csv_json)
        } catch (ValidationException e) {
            println ''
            log.error 'ERROR: Validation of input csv failed!'
            JSONObject exceptionJSON = e.toJSON()
            printExceptions(exceptionJSON, csv_json, log, enums)
            // log.error exceptionJSON.toString(2)
            println ''
            has_error = true
        }

        if (has_error) {
            System.exit(1)
        }
    }
    
    //
    // Loop over nested exceptions and print the causingException
    //
    private static void printExceptions(ex_json, params_json, log, enums, limit=5) {
        def causingExceptions = ex_json['causingExceptions']
        if (causingExceptions.length() == 0) {
            def m = ex_json['message'] =~ /required key \[([^\]]+)\] not found/
            // Missing required param
            if (m.matches()) {
                log.error "* Missing required column: ${m[0][1]}"
            }
            // Other base-level error
            else if (ex_json['pointerToViolation'] == '#') {
                log.error "* ${ex_json['message']}"
            }
            // Error with specific param
            else {
                def (_,row,column) = (ex_json['pointerToViolation'] =~ /^#\/(\d+)\/(.+)/)[0]
                int row_index = row.toInteger()
                int row_number = row_index + 2
                def param_val = params_json[row_index][column].toString()
                if (enums.containsKey(column)) {
                    def error_msg = "* Row ${row_number} column ${column}: '${param_val}' is not a valid choice (Available choices"
                    if (enums[column].size() > limit) {
                        log.error "${error_msg} (${limit} of ${enums[column].size()}): ${enums[column][0..limit-1].join(', ')}, ... )"
                    } else {
                        log.error "${error_msg}: ${enums[column].join(', ')})"
                    }
                } else {
                    log.error "* Row ${row_number} column ${column}: ${ex_json['message']} (${param_val})"
                }
            }
        }
        for (ex in causingExceptions) {
            printExceptions(ex, params_json, log, enums)
        }
    }
}