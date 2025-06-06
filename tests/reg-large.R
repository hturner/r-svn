#### Regression Tests that need "much" memory
#### (and / or  are slow even with enough GBytes of memory)

print(si <- sessionInfo(), locale=FALSE)
Sys.info()
## Run (currently _only_)  when inside tests/  by
'
time   make test-Large
' # giving ~ 35 min [R-devel 2019-01]

## From CRAN package 'sfsmisc':
Sys.memGB <- function (kind = "MemTotal")
{
    mm <- drop(read.dcf("/proc/meminfo", fields = kind))
    if (any(is.na(mm))) {
        warning("Non-existing 'kind': ", names(mm)[is.na(mm)][1])
        0
    } else if (!all(grepl(" kB$", mm)))  {
        warning("Memory info ", dQuote(kind),
                " is not returned in 'kB' aka kiloBytes")
        0
    } else
        as.numeric(sub(" kB$", "", mm))/(1000 * 1024)
}

availableGB <-
    if(file.exists("/proc/meminfo")) { # e.g. on Linux
	Sys.memGB("MemAvailable")
    } else {
	0 # unless we add something better here
    }
cat("Available (processor aka CPU) memory: ", round(availableGB, 1),
    "GB (Giga Bytes)\n")

if(.Machine$sizeof.pointer < 8) {
    cat(".Machine :\n"); str(.Machine)
    cat("not a 64-bit system -- forget about these tests!\n")
    q("no")
}

### Testing  readLines()  *large* file with embedded nul aka `\0'
##
## takes close to one minute and ~ 10 GB RAM
if(availableGB > 11) local(withAutoprint({
    ## File construction originally by Bill Dunlap, Cc: R-help,
    ##    Subject: Re: [R] readLines without skipNul=TRUE causes crash
    ##    Date: Mon, 17 Jul 2017 08:36:55 -0700
    tf <- tempfile(); file <- file(tf, "wb")
    txtLine <- c(rep(as.raw(32:127), 2^5), charToRaw("\n")) # <- want many lines
    system.time({
        for(i in 1:(2^15-1)) writeBin(rep_len(txtLine,    2^16), file)
        for(i in 1:(2^15-1)) writeBin(rep_len(as.raw(0L), 2^16), file)
    })
    close(file)
    log2(file.size(tf)) ## 31.99996
    ## now, this gave a segmentation fault, PR#17311 :
"FIXME: on 32-bit Linux (F 24), still see
Program received signal SIGSEGV, Segmentation fault.
... in do_readLines (call=0x8.., op=0x8.., ....)
    at ../../../R/src/main/connections.c:3852
3852		    if(c != '\n') buf[nbuf++] = (char) c; else break;
"
  if(.Machine$sizeof.pointer > 4) withAutoprint({
    system.time( x <- readLines(tf) ) # depending on disk,.. takes 15-50 seconds
    ##                ---------
    str(ncx <- nchar(x, "bytes"))
    ## int [1:688108] 3072 3072 3072 3072 3072 3072 3072 3072 ...
    tail(ncx) # ... 3072 3072 3072 1003
    table(ncx) # mostly 3072, then some 4075 and the last one
    head(iL <- which(ncx == 4075))
    stopifnot(diff(iL) == 21)
  }) else cat("32-bit: still seg.faulting - FIXME\n")
}))
## + 2 warnings


### Testing PR#17992  c() / unlist() name creation for large vectors
## Part 1
if(availableGB > 21) system.time({
    res <- c(a=raw(2), raw(2^31-1))
}) ## 36--44 sec elapsed (ada-16, ~ 120 GB available) after fix
## In R <= 3.4.1, took  51 sec elapsed, and gave Error .. :
##  attempt to set index 18446744071562067968/2147483649 in SET_STRING_ELT
##
if(FALSE) { # object.size() itself is taking a lot of time!
    os <- object.size(res)
} else {
    os <- structure(19327353184, class = "object_size")
    print(os, units = "GB") # 18
}
if(exists("res")) rm(res)
gc(reset = TRUE) # for the next step

### Testing PR#17992  c() / unlist() name creation for large vectors
## Part 2 (https://bugs.r-project.org/show_bug.cgi?id=17292#c4):
if(availableGB > 37) system.time({
    res <- c(a = list(rep(c(b=raw(1)), 2^31-2), raw(2)), recursive=TRUE)
})
## 437 sec elapsed (ada-16, ~ 120 GB available) after fix; then ada-20: 566 sec
## In R <= 3.4.1, took  475 sec  elapsed, and gave Error .. :
##    could not allocate memory (2048 Mb) in C function 'R_AllocStringBuffer'
## ((and that error msg is incorrect because of int overflow))
if(exists("res")) withAutoprint({
str(res) # is fast!
## Named raw [1:2147483648] 00 00 00 00 ...
## - attr(*, "names")= chr [1:2147483648] "a.b" "a.b" "a.b" "a.b" ...
gc() # back to ~ 18.4 GB
rm(res)
})
gc(reset = TRUE) # for the next step

## Large string's encodeString() -- PR#15885
if(availableGB > 4) system.time(local(withAutoprint({
    txt <- strrep("test me:", 53687091); object.size(txt) # 429'496'824 bytes
    nc <- nchar(txt) ## NB this is larger than maximal integer:
    nc*5L+8L # NA + Warning   'NAs produced by integer overflow'
    en <- encodeString(txt)
    ## encodeString() seg.faulted in R <= 3.4.1
    stopifnot(identical(txt,en)) # encoding did not change simple ASCII
})))
## 52 sec elapsed [nb-mm4, 8 GB]; then 66.7 [ada-20; much more GB]


## pretty(x, n) for n = <large> or  large diff(range(x) gave overflow in C code
if(availableGB > 6) system.time(withAutoprint({
    r <- pretty(c(-1,1)*1e300, n = 449423288, min.n = 1)
    head(r) ; length(r) # was only 21 in  R < 3.5.0
    stopifnot(all.equal(length(r), 400000001, tol = 0.1))
})) ## 4.8--5.5 sec.
rm(r)
gc()

n <- 4e4 # << for quick testing, comment next line
n <- 2.2e9

if(availableGB > 60) withAutoprint({
    n/.Machine$integer.max  # 1.024 ==> need  long vectors!
    ii <- seq_len(n)          #   user  system elapsed  [seq_len() fast: ALTREP "compact"]
    system.time(ii <- ii + 0) #  6.726  17.558  24.450 (slow!, seen faster)
    system.time(i2 <- ii[-n]) # 14.267  23.532  37.918 (slow!, seen slower: el.= 51)
    ##
    ## NB: keep n, ii, i2 for "below"
})
## In R <= 3.4.1 :
## Program received signal SIGSEGV, Segmentation fault.
## 0x00000000005a0daf in realSubscript (call=0x3f01408, stretch=<optimized out>,
##     nx=2200000000, ns=1, s=0x426db18) at ../../../R/src/main/subscript.c:691
## 691			    LOGICAL(indx)[ix] = 0;

if(availableGB > 99) withAutoprint({
    system.time( x <- ii/n )            #   5.45 user; 11.5--14.36 elapsed
    system.time( y <- sin(pi*x) )       #  42 user; 48.9--..  elapsed
    system.time(sorted <- !is.unsorted(x)) # ~ 4 elapsed
    stopifnot(sorted)
    ## default n (= "nout") = 50:
    system.time(ap1 <- approx(x,y, ties = "ordered"))# 15 user; 25 elapsed
    stopifnot(exprs = {
	is.list(ap1)
	names(ap1) == c("x","y")
	length(ap1$x) == 50
	all.equal(ap1$y, sin(pi*ap1$x), tol= 1e-9)
    })
    rm(ap1) # keep x,y,n,ii,i2
    gc()     # --> max used: 92322 Mb
})

## which() and ifelse() working for long vectors
if(availableGB > 165) withAutoprint({
    system.time(iis <- which(isMl <- ii < 9999)) # 5.8 user,  8.8 elapsed
    gc() # 59 GB max used
    system.time(r <- ifelse(isMl, ii, ii*1.125)) #        user  system elapsed
    stopifnot(exprs = {                 # in R 3.5.2 : 124.989 174.726 300.656
	## GB's ifelse() + using which(<long>) 3.6.0 :  71.815  81.823 154.124
	length(r) == n
        iis == seq_len(9998)
    })
    rm(isMl, iis, r)
})
gc() # 159 GB max used

if(availableGB > 211) withAutoprint({ ## continuing from above
    ## both large (x,y) *and* large output (x,y):
    system.time(xo <- x + 1/(2*n))     # ~ 9 elapsed
    system.time(ap <- approx(x,y, ties = "ordered", xout = xo))
                                       # 194 user, 214--500 elapsed
    gc(reset = TRUE) # showing max.used ~ 1..... Mb
    stopifnot(exprs = {
	is.list(ap)
	names(ap) == c("x","y")
	length(ap$x) == n
	is.na(ap$y[n]) # because ap$x[n] > 1, i.e., outside of [0,1]
	all.equal(ap$y[i2], sin(pi*xo[i2]), tol= if(n < 1e7) 1e-8 else 1e-15)
    })
    rm(ap); gc() # showing used 83930 Mb | max.used 210356.6 Mb
    ## only large x,y :
    system.time(apf <- approxfun(x,y, ties="ordered", rule = 2))# elapsed: ~26s
    xi <- seq(0, 1, by = 2^-12) ## linear interpol. is less accurate than spline:
    stopifnot(all.equal(apf(xi), sin(pi*xi), tol= if(n < 1e7) 1e-7 else 1e-11))
    rm(apf); gc() # (~ unchanged)
    system.time(ssf <- splinefun(x,y, ties = "ordered"))
                                        # elapsed 120 s; using ~ 158 GB
    system.time(ss  <- spline   (x,y, ties = "ordered", xout = xi))
                                        # elapsed 126--265 s; using ~ 207 GB
    gc()
    stopifnot(exprs = {
	is.list(ss)
	names(ss) == c("x","y")
	length(ss$y) == length(xi)
	all.equal(ss$y   , sin(pi*xi), tol= 1e-15)
	all.equal(ssf(xi), ss$y,       tol= 1e-15)
    })
    rm(x, y, xo, ss, ssf) # remove long vector objects
    gc(reset=TRUE)
})

## sum(<Integer|Logical>) -- should no longer overflow: ----------------------------------------
## 1) sum(<long logical>) == counting
if(availableGB > 24) withAutoprint({
    system.time(L <- rep.int((0:15) %% 7 == 2, 2^28))# -> length 2^32; ~ 22 sec
    print(object.size(L), unit="GB") # 16 GB
    system.time(sL <- sum(L)) # 8.4 sec
    stopifnot(exprs = {
        is.logical(L)
        length(L) == 2^32
        !is.integer(length(L))
        is.integer(sL)
        identical(sL, as.integer(2^29))
    })
}) ## sL would be NA with an "integer overflow" warning in R <= 3.4.x
gc(reset=TRUE)

## 2) many (and relatively long and large) integers
L <- as.integer(2^31 - 1)## = 2147483647L = .Machine$integer.max ("everywhere")
## a "small" example with this is in ./reg-tests-1d.R (see 'x24')
if(availableGB > 12) withAutoprint({
    system.time(x31 <- rep.int(L, 2^31+1)) # sum = 2^62 - 1 =.= 2^62 // ~ 5.5 sec
    print(object.size(x31), unit = "GB") # 8 G
    system.time(S <- sum(x31)) # ~ 2 sec
    system.time(S.4 <- sum(x31, x31, x31, x31)) # 8 sec
    stopifnot(is.integer(x31),
              identical(S,   2^62),
              identical(S.4, 2^64))
    system.time(x32 <- c(x31, x31)) # 13 user | 20.8 elapsed  (and 16 GB)
    rm(x31)# now,  sum vvv  will switch to use irsum() [double accumulator]
    system.time(S.2 <- sum(x32)) # 8 sec
    stopifnot(S.2 == 2^63)
    rm(x32)
})


## seq() remaining integer: (PR 17497, comment #9)
if(availableGB > 16) withAutoprint({
    i <- as.integer(2^30)
    system.time(i2.31 <- seq(-i, by=1L, length=2*i+1)) # 11.1 user | 19.2 elapsed
    object.size(i2.31) # 8'589'934'648 bytes [ was 17.17 GB in R <= 3.5.x ]
    stopifnot(is.integer(i2.31),  i2.31[1] == -i,  i2.31[length(i2.31)] == i)

    ## pmax(), pmin() with long vectors, PR 17533
    if(availableGB > 24) withAutoprint({
        system.time(i2.31 <- pmin(i2.31, 0L)) # 7.2 sec user | 11.2 elapsed
        str(i2.31)
        system.time(stopifnot(i2.31[(i+1):length(i2.31)] == 0)) # 16.7 user | 28.0 elapsed
    })
})


## match(<long character>, *)  PR#17552
if(availableGB > 44) withAutoprint({ ## seen 40 G ('RES')
    system.time(m <- match(rep("a", 2^31), "a")) # 34.7 sec user (55 elapsed)
    stopifnot(all(m == 1L))
    rm(m)
    system.time({x <- character(2^31); x[26:1] <- letters }) # 1.6 user | 9.4 elapsed
    system.time(m <- match(x, "a"))# 18.2 user | 51.6 elapsed
    head(m, 30)
    system.time(stopifnot(m[26] == 1L, is.na(m[-26])))
    rm(x, m)
})


## readBin() and writeBin() for long rawConnection s, PR#17665
## -------       --------            -------------
if(availableGB > 14) withAutoprint({ ## seen 11.6 G
    vec <- rep(0, 3e8) # object.size(vec) > 2^31
    raw_con <- rawConnection(serialize(vec, NULL)) # ~ 5 sec.
    ## Stepping through this connection gives an error after the 2^31st element:
    repeat {
        x <- readBin(raw_con, "raw", n = 1e+06)
        if(length(x) == 0)
            break
        cat(".")
    }; cat("\n")
    ## Error in readBin(raw_con, "raw", n = 1e+06) : too large a block specified
})

## writeBin() for long vectors
if(availableGB > 20) withAutoprint({ ## seen 20.9 G
    x <- raw(2^31)
    writeBin(x, con = nullfile())

    con <- rawConnection(raw(0L), "w")
    writeBin(x, con = con)
    stopifnot(identical(x, rawConnectionValue(con)))

    system.time(x <- pi*seq_len(2.1*2^30)) # 25 sec
    zzfil <- tempfile("test-large-bin")
    zz <- file(zzfil, "wb") ## file size will be 2.5 GB !!!
    system.time(z <- writeBin(x, zz)) # 32 sec
    stopifnot(is.null(z))
    close(zz); zz <- file(zzfil, "rb")
    system.time(r <- readBin(zz, double(), n = length(x) + 999)) # 32 sec
    system.time(stopifnot(identical(x, r))) # 24 sec
    close(zz); rm(r, zz)
})


## predict(loess(.), se=TRUE) for "large" sample size -- PR#17121
## No need for very much memory, but is slow and should do several ex.
mkDat <- function(n) {
    x <- 5*(1:n)/(n+1)
    data.frame(x = x, y = sin(pi*x^2) * exp(-x/2) + rnorm(n)/8)
}
set.seed(1); dat <- mkDat(n = 42000)
system.time( # 12.7 sec (lynne ~ 2021)
    fit <- loess(y~x, data=dat)
)
r <- tools::assertError(
   predict(fit, newdata=data.frame(x=.5), se=TRUE)
 , verbose=TRUE) #
## typically would not seg.fault but give Calloc(..) error (with *wrong* size)
stopifnot(grepl("^workspace .* is too large .* 'se = TRUE'", r[[1]]$message))


## PR#17330 :  '[[<-' for index 2^31 :
(i <- 2^31) > .Machine$integer.max
system.time(x <- raw(i)) # ~ 0.8 sec ; needs 2 GB
x [i]  <- r1 <- as.raw(1); stopifnot(x [i]  == r1)
x[[i]] <- r2 <- as.raw(2); stopifnot(x[[i]] == r2)
x[[i]] <- r3 <- as.raw(3); stopifnot(x[[i]] == r3)
## last two failed in R <= 4.0.n {even with large vectors}


## print()ing {up to max.print only!} of long vectors;
## including named and "generic" (= list):
stopifnot((n <- 2^31 + 352) > .Machine$integer.max)
system.time(L <- integer(n))        #   5.8 sec {ada-20}
system.time(LL <- vector("list", n))# ~15   sec {ada-20}
system.time(nm <- c(LETTERS, letters, rep("xx", length(L) - 2*26)))
## between 55 and 76 sec {ada-20, 2022-01-07}  user  system elapsed
Ln <- L
## FIXME? takes about 2 secs, but these are *not* seen by system.time (!!)
system.time(names(Ln) <- nm)
## user  system elapsed
##    0       0       0
op <- options(max.print = 300)
L
## now (after using %lld) gives
## [ reached getOption("max.print") -- omitted 2147483700 entries ]
## before, it gave   ..... -- omitted -2147483596 entries
##                                   ^^^
Ln # gave  Error: long vectors not supported yet: ...
LL # gave  Error: long vect...
options(op)
rm(Ln, L)

## PR#17977 --- x[<fractional>] behavior should fulfill x[i] === x[as.integer(i)]
## large (no overflow in index computations!) -- needs 2 GB
LL <- matrix(as.raw(1:2), 2, 2^30)
ca.half <- 0.5+ (eps <- unique(sort(outer(2^-c(16, 21, 26, 30), -1:1))))
print(eps, digits=3)
LL[cbind(2, ca.half)]   # should be of length 0, too: ca.half ~= 0.5
LL[cbind(1, 1+ca.half)] # should be constantly == raw(1L) '01'
LL[cbind(2+ca.half, 1)] # all 02 -- failed in R <= 4.1.x
LL[cbind(-ca.half, 1)]  # raw(0) --  "      "    "
stopifnot(exprs = {
    length(LL[cbind(2, ca.half)]) == 0
    LL[cbind(1, 1+ca.half)] == as.raw(1L)
    LL[cbind(2+ca.half, 1)] == as.raw(2L)
    length(LL[cbind( -ca.half, 1)]) == 0
})

if(availableGB > 10) withAutoprint({ ## PR#18612
    ##  Summary: integer overflow in matrix(<long vector>, nrow, ncol)
    ##           due to wrong format specifiers
    ## Reporter: Mikael Jagan @ McMaster
    M <- .Machine$integer.max
    x <- raw(M + 1)
    y <- raw(2 * M)
    (m1 <- conditionMessage(tryCatch(matrix(x, M, 1L), warning = identity)))
    ## was "data length [-2147483648] is not a sub-multiple or multiple ...."
    (m2 <- conditionMessage(tryCatch(matrix(x, 1L, M), warning = identity)))
    ## was "data length [-2147483648] is not a sub-multiple or multiple ...."
    (m3 <- conditionMessage(tryCatch(matrix(y, M, 1L), warning = identity)))
    ## was "data length differs from size of matrix: [-2 != 2147483647 x 1]"
    (m4 <- conditionMessage(tryCatch(matrix(y, 1L, M), warning = identity)))
    ## "data length differs from size of matrix: [-2 != 1 x 2147483647]"
    ## triggering those in dimsgets():
    (m5 <- conditionMessage(tryCatch(dim(y) <- c(1L, M), error = identity)))
    (m6 <- conditionMessage(tryCatch(dim(x) <- c(M, 1L), error = identity)))
    stopifnot(exprs = {
        grepl(paste0("data length [", M+1, "] is not"),   c(m1,m2), fixed=TRUE)
        grepl(sprintf("size of matrix: [%.0f != ", 2*M),  c(m3,m4), fixed=TRUE)
        grepl(paste0("dims [product ", M, "] do not "),   c(m5,m6), fixed=TRUE)
    })
    rm(y)
})

x <- -1:2^31 # (immediate: ALTREP)
system.time( r <- rank(x) )# Error about invalid length() -- PR#18617, in R <= 4.3.2
## seen 260 sec (!)
head(r)
stopifnot(r[1:6] == 1:6)
rm(r,x)

## rank(<largish) -- PR#18630
vals <- 1:1475000000 # (immediate, thanks to ALTREP)
system.time( ranks <- rank(vals) ) # ~ 130 sec
head(ranks, 11)
stopifnot(ranks[1:11] == 1:11, min(ranks) == 1)
## min(ranks) was -1073741824, in R <= 4.3.2



gc() # NB the "max used"
proc.time() # total  [ ~ 50 minutes in full case, 2024-01; was 40' in 2019-04-12]
