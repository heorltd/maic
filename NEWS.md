News for package maic
================
Rob Young (<rob.young@heor.co.uk>)
17/04/2020

# 11 May 2020

The following fixes have been made:

-   Reporting when probabilities are equal to 1 or 0 no longer fails
    (\#1)
-    is now used to summarise
-   Index and dictionary inputs are forced to data.frame (\#2)
-   Clarification that the p-values for difference in mean are
    stochastically derived (\#3)

# 17 April 2020

This is the first CRAN upload of . Read the documentation, check the
function signatures, have a look at the attributes of what gets
returned. There are many little features that have been added to improve
your workflow, but there are major changes from the development
versions, especially:

-   Weights will now always be same same length as the number of rows in
    your data. You will no longer need to manually put in zeros for
    those subjects excluded due to lack of overlap.
