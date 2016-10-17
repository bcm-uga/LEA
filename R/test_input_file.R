test_input_file <- function(name, extension)
{
    # list of possible extension
    l = c("geno", "lfmm", "vcf", "ancestrymap", "ped")

    # obtain the extension of name
    ext = getExtension(basename(name))

    # init input_file
    input_file = name;
    # check if the exstension is known
    if (!(ext %in% l)) {
        p = paste("The extension (.", ext,") is unknown (file \"",name,"\").\n",
            "Please, use one of the following format extension: .geno,", 
            " .lfmm, .vcf, .ancestrymap, .ped.", sep="")
        stop(p);
    # if not the correct format, convert
    } else if (ext != extension) {
        input_file = setExtension(name, paste(".", extension, sep=""))
        print("*********************************************************");
        print(paste(" Conversion from the ", ext," format to the ", extension,
            " format", sep = ""));
        print("*********************************************************");
    
        if (extension == "lfmm") {    
            if(ext == "geno") {
                input_file = geno2lfmm(name, force = FALSE)
            } else if (ext == "ancestrymap") {
                input_file = ancestrymap2lfmm(name, force = FALSE)
            } else if (ext == "vcf") {
                input_file = vcf2lfmm(name, force = FALSE)
            } else if (ext == "ped") {
                input_file = ped2lfmm(name, force = FALSE)
            }
        } else if (extension == "geno") {
            if(ext == "lfmm") {
                input_file = lfmm2geno(name, force = FALSE)
            } else if (ext == "ancestrymap") {
                input_file = ancestrymap2geno(name, force = FALSE)
            } else if (ext == "vcf") {
                input_file = vcf2geno(name, force = FALSE)
            } else if (ext == "ped") {
                input_file = ped2geno(name, force = FALSE)
            }
        }
    }

    return(input_file);
}
