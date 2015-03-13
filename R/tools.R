
start.silent <- function (silent, d, extension)
{
    if (silent) { 
        log = paste(d, extension, ".log", sep="")
        sink(log)
    }
}

stop.silent <- function(silent)
{
    if (silent)
        sink()
}

modified <- function(file1, file2)
{
    if (! file.exists(file2) || !file.exists(file1)) {
        return(TRUE)        
    } else {
        return (file.info(file1)$mtime > file.info(file2)$mtime)
    }
}

projectLfmmLoad <- function(input.file, environment.file, project)
{
    inp = basename(setExtension(normalizePath(input.file), ""))
    var = basename(setExtension(normalizePath(environment.file), ""))
    projectName = paste(dirname(normalizePath(input.file)), "/",
        inp, "_", var, ".lfmmProject", sep="")
    # load the project
    if (!project == "new" && file.exists(projectName)) {
        proj = load.lfmmProject(projectName)
        if (proj@creationTime < file.info(input.file)$mtime) {
            stop("The input file has been modified since the creation of the",
                " project.\nIf the input file is different, the results ",
                "concatenating all runs can be false.\nTo remove the current",
                " project and start a new one, add the option 'project = ",
                "new'.\nTo continue with the same project, add the option ",
                "'project = force'.")
        }
        if (proj@environment.file != environment.file) {
            stop(paste("A project with the same variable ('",var,"') and a",
                "different path ('", proj@environment.file, "'and '", 
                environment.file, " already exists. Please, set a different",
                " name for the variable file.", sep=""))
        }
        if (proj@creationTime < file.info(environment.file)$mtime) {
            stop("The variable file has been modified since the creation of ",
                "the project.\nIf the variable file is different, the results",
                " concatenating all runs can be false.\nTo remove the current",
                " project and start a new one, add the option 'project = ", 
                "new'.\nTo continue with the same project, add the option ",
                "'project = force'. To create a new project and keep the ",
                "previous one, change the name of the variable file.")
        }
    # create a new project
    } else {
        if (file.exists(projectName)) {
            remove.lfmmProject(projectName);
        }
        proj = new("lfmmProject")
        proj@creationTime = Sys.time()
        #files 
        proj@input.file = normalizePath(input.file)
        proj@environment.file = normalizePath(environment.file)
        # directory
        proj@directory = setExtension(projectName, ".lfmm/")
        unlink(proj@directory, recursive = TRUE)
        dir.create(proj@directory, showWarnings = FALSE)
        # lfmmProject Files
        proj@lfmmProject.file = projectName
        save.lfmmProject(proj)

    }
    return(proj)
}

projectSnmfLoad <- function(input.file, project)
{

    projectName = setExtension(paste(dirname(normalizePath(input.file)), "/", 
        basename(input.file), sep=""), ".snmfProject")
    # load the project
    if (!project == "new" && file.exists(projectName)) {
        proj = load.snmfProject(projectName)
        # file input file not modified since the start of the project
        if (proj@creationTime >= file.info(input.file)$mtime 
            || project == "force") {
            return(proj)
        } else {
            stop("The input file has been modified since the creation of the",
                " project.\nIf the input file is different, the results ",
                "concatenating all runs can be false.\nTo remove the current",
                " project and start a new one, add the option 'project = ",
                "new'.\nTo continue with the same project, add the option ",
                "'project = force'.")
        }
    # create a new project
    } else {        
        if (file.exists(projectName)) {
          remove.snmfProject(projectName);
        }
        proj = new("snmfProject")
        proj@creationTime = Sys.time()
        # files
        proj@input.file = normalizePath(input.file)
        # directory    
        proj@directory = setExtension(paste(dirname(normalizePath(input.file)), 
            "/", basename(input.file), sep=""), ".snmf/")
        unlink(proj@directory, recursive = TRUE)
        dir.create(proj@directory, showWarnings = FALSE)
        # snmfProject.file
        proj@snmfProject.file = projectName
        # masked
        dir.create(paste(proj@directory, "masked/", sep=""), 
            showWarnings = FALSE) 
        save.snmfProject(proj)

        return(proj)
    }
}

getExtension <- function(file)
{
    l = strsplit(file, "\\.")[[1]]
    return(l[length(l)])
}

setExtension <- function(file, ext)
{
    out = dirname(file)
    l = strsplit(basename(file), "\\.")[[1]]
    if (length(l) >= 2) {
        out = paste(out, l[1], sep="/")
        if (length(l) >= 3) {
            for (i in 2:(length(l)-1)) 
                out = paste(out, l[i], sep=".")
            }
            out = paste(out, ext, sep="")
        } else if (length(l) == 1) {
            out = paste(out, l[1], sep="/")
            out = paste(out, ext, sep="")
        } else {
        out = NULL
    }

    return(out)     
}

currentDir <- function(file)
{
    substr(file, nchar(getwd())+2, nchar(file))
}
