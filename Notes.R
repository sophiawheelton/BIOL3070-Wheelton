# Day 1 Example Plot

#import the file:
viremia <- read.csv("viremia_data_full.csv")

#view the data in a note:
View(viremia)

#name columns
colnames(viremia) <- c("Bird","n","Family","Order","1","3","4","6")

#choose some colors
cols <-c("black","grey",rainbow(26)[4:26])

#Plot by species - did not work
plot(c(1,3,4,6),as.numberic(viremia[1,6:9]),
  type = "l",lwd =2, ylim = range(viremia [,6:9],na.rm=TRUE),
  xlab = "Day Postinfection",
  ylab = "Log PFU/ml Serum")
for (i in 2:nrow(viremia)){
  lines(c(1,3,4,6),as.numeric(viremia[i,6:9]),lwd =2, col = cols)
}


------------------------------------------------------------------
# Day 2 Objects and Data Types

# an object stores information in R (an object oriented language)
x <- 5 # numeric
y <- "cat" # character - object y is assigned a character value "cat"
z <- TRUE # logical - only TRUE (T) or FALSE (F) 

# code is always from top to bottom, left to right

# numeric = numbers
# character = text
# logical = TRUE/FALSE 
# factor = fixed categories (ie sex male/female, pets cat/dog/fish/bir)

#data types change how R interprets and uses data (ie "5"+1->error - character can't add to number)

#console: where code runs immediately
#script (.Rmd or .R): where code is written, saved, reused
#best practice: test in console -> save in script/RmD

#Practice:
> a<-5
> b<-"5"
> c<-TRUE
> d<-factor(c("A","B","A","C"))
> a+2
[1] 7
> b == 5 #just checking what was already stated/if a given condition is true
[1] TRUE
> !c #Bang operator (!) means everything but or the opposite
[1] FALSE
> levels(d)
[1] "A" "B" "C"

#"I use AI to help with error correction" - how we are allowed to use AI 
#can't put spaces in object names

# is called commenting

----------------------------------------------------------------------  
#Day 3

#Rmd stands for R mark down
#reproducible reports=text+code+output - a file that saves text
#headers (# in markdown)
#code chunks ('''{r})
#commenting in/out of code

#write code in grey and text outside
# "##" is headers
#outline shows you headers
#Knit "prints" it
  
#output in the top section changes how it is downloaded (PDF, word, HTML)

---------------------------------------------------------------------------

#Day 4
# <- and = are the same, == checks and results in "true"
#for factor, you have to say factor()
# "4" is different than 4 which is different than 4.0 - testing them against each other would result in FALSE

#R Markdown
  #a format for writing reproducible, dynamic reports with R - use it to embed code and results into slideshows, pdfs, etc
  #open, write, embed and render (replace R code with its output and transform the report into whatever you want)
  #menu: File - New File - R Markdown - class of output - OK
  #YAML header is title, author, output, etc.
  #use # for headers: # is big then ## smaller then ### smallest
  #Ctrl+Alt+I puts in 3 ticks (code chunks), subract to 1 for incline code
    #code chunks show the code in the final report
    #inline code replaces the code with the results in the final report
  #commenting outside of code chunk only adds the #, if you want it to be emitted in final report: <!-- your comment -->
  #in the code chunk, the comment will be visible with the #
  #Ctrl+Shift+C is shortcut for commenting
  #eval(or echo)=FALSE - this will leave out the code but will show the output

  
#Day 5
  #syntax
  my_function<-function(parameters) {
    #code to execute.
  }
  my_function(parameters)
  
  odd_or_even <- function(num) {
  }
    
  ?round()
  ?c()
  ?factor()
  ?plot()
  
  #? is help menu, made by random people
  #round
  #x and digits (arguments) are required - goes in ()
  #default is what it decides you want if you don't specify with an argument
  round(pi) #or
  round(x=3.141593) #define x and leave digits default (0)
  round(x=3.141593,digits=3) #digits are after the decimal points
  #homework: go through helps for each of the above and understand them
  
  
  #Day 6 
  #practice for Day 5
  round(4,digits=3)
  round(4.39872, digits=3)
  ?round()
  round(4.39872)
  ?mean  
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,NA,2))
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,2))
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,2),trim=0.1)
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,NA,2),na.rm=TRUE)

#Day 7
  #read.table(), read.csv() -> import data
  #write.table(), write. csv() -> save/export
  #file paths (relative vs absolute)
  #tip: always check your data after import
    #head() peek at first rows
    #view() open spreadsheet view
    #str() 
  
  
?read.csv #typically use commas as separaters
data <- read.csv(file="viremia_data_full.csv")

head(data) #headers
View(data) #needs to be capitalized
str(data) #structure
dim(data) #column and row dimensions

data <- read.csv(file="bloodmeal_for_BIOL3070.csv")

#Day 8 Mini-report Rubric
#total: 25 points - five grading categories, 5 points each
#1. Structure and Completion
  #All sections present (Abstract, Background, Question/Hypothesis, Methods, Discussion, Conclusion, References)
  #Abstract presents and covers full project (written last)
  #Knits without errors
#2. Background and Question Framing
  #Clear background with at least 1 cited scientific reference
  #Connects to course material (e.g. viremia plot)
  #States study question, hypothesis, and prediction clearly
#3. Methods and Analysis
  #Methods written in prose, not just code
  #At lease 2 analyses/plots with reproducible R code
  #Code commented and organized
#4. Interpretation and Discussion 
  #Interpret each analysis/plot in text
  #Tie results back to hypothesis
  #Note at least 1 limitation or uncertainty
#5. Clarity, References, and Reproducibility
  #Clear, professional writing and formatting 
  #All references cited (scientific + AI if used)
  #Knitted .md file shared via GitHub URL
  #Report looks polished and reproducible

#Day 9 - Various
#commit saves a version of the work
  #records what changed along with a short message
#push uploads the save changes (commits) from RStudio/Posit Cloud to your GitHub repository online
  #how you publish your newest version
