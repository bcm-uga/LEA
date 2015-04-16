create.dataset <- function(input.file, output.file, seed = -1, 
    percentage = 0.05) 
{
    # test arguments and init
    # input file
    input.file = test_character("input.file", input.file, NULL)
    # check extension and convert if necessary
    input.file = test_input_file(input.file, "geno")
    # output file    
    tmp = sub("([^.]+)\\.[[:alnum:]]+$", "\\1_I.geno",input.file)
    output.file = test_character("output.file", output.file, tmp)
    # seed
    seed = test_integer("seed", seed, 
        as.integer(runif(1)*.Machine$integer.max))
    if (seed == -1) 
        seed = as.integer(runif(1)*.Machine$integer.max)
    set.seed(seed) # init seed
    print(seed)
    # percentage
    percentage = test_double("percentage", percentage, 0.05)
    if (percentage <= 0 || percentage >= 1)
        stop(paste("'percentage' argument has to be of type double and ", 
            "between 0 and 1.", sep=""))

    print("*************************************");
    print("*          create.dataset            *");
    print("*************************************");

    # run method
    .C("R_createDataSet", 
    as.character(input.file),
    as.integer(seed),
    as.double(percentage),
    as.character(output.file)
    );

    # create output 
    return(output.file);
}
