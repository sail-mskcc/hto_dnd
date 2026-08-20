version 1.0

task Exception {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Raise an exception with a custom message"
    }

    parameter_meta {
        message: { help: "The exception message to raise." }
    }

    input {
        String message
    }

    command <<<
        echo "~{message}" >&2
        exit 1
    >>>

    output {
        String message_out = message
    }

    runtime {
    }
}