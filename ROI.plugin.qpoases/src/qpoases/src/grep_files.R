

files <- list.files(pattern = ".cpp")

x <- unlist(lapply(files, readLines))
i <- grep("return ", x)

y <- x[i]
m <- gregexpr("[A-Z_]*", y)

z <- unlist(regmatches(y, m))
z <- z[nchar(z) > 0]

tab <- table(z)

tab[nchar(names(tab)) <= 10]
length(tab[nchar(names(tab)) > 10])

grep("success", names(tab), ignore.case = TRUE, value = TRUE)



