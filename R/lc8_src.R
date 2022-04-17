## LC8 package source
## source: https://github.com/LiChenPU/Formula_manipulation/
## Reference: https://www.nature.com/articles/s41592-021-01303-3
##



my_calculate_formula = function(Formula1,Formula2,sign = 1, Valid_formula = FALSE){
  formula2_ls = list()
  for(i in 1:length(Formula2)){
    formula2 = Formula2[i]
    {
      formula2 <- gsub("D", "[2]H", formula2)
      ende2 <- nchar(formula2)
      element2 <- c()
      number2 <- c()
      j <- c(1)
      while (j <= ende2) {
        if (substr(formula2, j, j) == c("[")) {
          b <- j
          while (any(substr(formula2, j, j) == c("]")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element2 <- c(element2, substr(formula2, b, m))
        }
        if (any(substr(formula2, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) != TRUE) {
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element2 <- c(element2, substr(formula2, k, m))
        }
        if (any(substr(formula2, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) == TRUE) {
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) ==
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number2 <- c(number2, as.numeric(substr(formula2,
                                                  k, m)))
        }
        j <- j + 1

      }

    }
    number2 = number2 * sign
    formula2_ls[[length(formula2_ls)+1]] = list(element2, number2)
  }

  formula1_ls = list()
  for(i in 1:length(Formula1)){
    formula1 = Formula1[i]
    {
      formula1 <- gsub("D", "[2]H", formula1)
      ende2 <- nchar(formula1)
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= ende2) {
        if (substr(formula1, j, j) == c("[")) {
          b <- j
          while (any(substr(formula1, j, j) == c("]")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element1 <- c(element1, substr(formula1, b, m))
        }
        if (any(substr(formula1, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) != TRUE) {
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formula1, k, m))
        }
        if (any(substr(formula1, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) == TRUE) {
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) ==
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formula1,
                                                  k, m)))
        }
        j <- j + 1
      }
    }
    formula1_ls[[length(formula1_ls)+1]] = list(element1, number1)
  }

  formula_mat = matrix("", nrow = length(formula1_ls), ncol = length(formula2_ls))
  valid_mat = matrix(F, nrow = length(formula1_ls), ncol = length(formula2_ls))
  for(row_i in 1:length(Formula1)){
    element1 = formula1_ls[[row_i]][[1]]
    number1 = formula1_ls[[row_i]][[2]]
    for(col_j in 1:length(Formula2)){
      element2 = formula2_ls[[col_j]][[1]]
      number2 = formula2_ls[[col_j]][[2]]
      count_all = number1
      elem_all = element1
      for(i in 1:length(element2)){
        if(any(element2[i]==elem_all) == TRUE){
          count_all[elem_all==element2[i]] = count_all[elem_all==element2[i]]+number2[i]
        }else{
          count_all = c(count_all,number2[i])
          elem_all = c(elem_all,element2[i])
        }
      }

      if(any(count_all<0)){
        valid_mat[row_i,col_j]=F
      }

      elem_order = order(elem_all)
      count_all = count_all[elem_order]
      elem_all = elem_all[elem_order]

      elem_all=elem_all[count_all!=0]
      count_all=count_all[count_all!=0]
      formula_all=c()
      for (i in 1:length(count_all)) {
        formula_all <- c(formula_all, elem_all[i],
                         count_all[i])
      }
      formula_mat[row_i,col_j]=paste0(formula_all, collapse = "")
    }
  }

  if(length(Formula1)==1 & length(Formula2)==1){
    formula_mat = as.vector(formula_mat)
    valid_mat = as.vector(valid_mat)
  }

  if(Valid_formula){
    return(list(formula_mat, valid_mat))
  } else {
    return(formula_mat)
  }
}



formula_mz <- function(formula = "C2H4O1", charge = 0,elem_table = lc8::elem_table){

  mz = numeric()
  for(i in 1:length(formula)){
    temp_formula = formula[i]
    if(is.na(temp_formula) | is.null(temp_formula)){
      mz[i] = NA
      next
    }
    if(temp_formula == ""){
      mz[i] = 0
      next
    }

    temp_formula <- gsub("D", "[2]H", temp_formula)
    ende2 <- nchar(temp_formula)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(temp_formula, j, j) == c("[")) {
        b <- j
        while (any(substr(temp_formula, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, b, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".","0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, k, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".","0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(temp_formula,
                                                k, m)))
      }
      j <- j + 1

    }

    mz[i]=0

    for(j in 1:length(element2)){
      mz[i] = mz[i] + elem_table$mass[element2[j] == elem_table$element] * number2[j]
    }


    e_mass = 0.00054857990943
    mz[i] = mz[i] - e_mass*charge

  }
  return(mz)


}


#' formula_rdbe
#'
#' @param formula e.g. "C2H4O1"
#' @param elem_table a table records unsaturation
#' @noRd
#' @return the ring and double bond number
#'
#' @examples formula_rdbe(formula = "C2H4O1", elem_table = lc8::elem_table)
formula_rdbe = function(formula = "C2H4O1", elem_table = lc8::elem_table){

  rdbe = numeric()
  for(i in 1:length(formula)){

    temp_formula = formula[i]
    temp_formula <- gsub("D", "[2]H", temp_formula)
    ende2 <- nchar(temp_formula)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(temp_formula, j, j) == c("[")) {
        b <- j
        while (any(substr(temp_formula, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, b, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, k, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(temp_formula,
                                                k, m)))
      }
      j <- j + 1

    }

    rdbe[i]=1
    for (j in 1:length(element2)) {
      rdbe[i] = rdbe[i] + elem_table$unsaturation[element2[j] ==
                                                    elem_table$element] * number2[j]
    }
  }
  return(rdbe)
}


#' isotopic_abundance
#'
#' @param formula e.g. "C2H4O1"
#' @param elem_table a table records unsaturation
#' @param isotope e.g. [13]C, or [13]C2 for M+2
#' @noRd
#' @return the ratio of given formula and its isotopic peak
#'
#' @examples isotopic_abundance(formula = "C2H4O1", elem_table = lc8::elem_table)
isotopic_abundance = function(formula = "C2H4O1", isotope = "[13]C1", elem_table = lc8::elem_table){
  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {

    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }


  number1 = number2
  element1 = element2


  formula <- gsub("D", "[2]H", isotope)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {

    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
        if(j>ende2){   break      }
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if(j>ende2){
      number2 = c(number2, 1)
      break
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
        if(j>ende2){   break      }
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if(j>ende2){
      number2 = c(number2, 1)
      break
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }

  selected_isotopes = grepl("\\[\\d+\\]",element2)
  number_isotopes = number2[selected_isotopes]
  element_isotopes = element2[selected_isotopes]

  ratio = 1
  for(i in 1:length(element_isotopes)){
    element_isotope = element_isotopes[i]
    element_parent = gsub("\\[\\d+\\]", "", element_isotope)
    num_isotope_formula2 = number_isotopes[i]


    num_parent_formula1 = 0
    num_isotpe_formula1 = 0
    if(any(element1 == element_parent)){
      num_parent_formula1 = number1[element1 == element_parent]
    }
    if(any(element1 == element_isotope)){
      num_isotpe_formula1 = number1[element1 == element_isotope]
    }

    natural_abundance = elem_table$abundance[elem_table$element == element_isotope]

    inten_formula1 = choose((num_parent_formula1+num_isotpe_formula1),num_isotpe_formula1)*natural_abundance^(num_isotpe_formula1)
    inten_formula2 = choose((num_parent_formula1+num_isotpe_formula1),(num_isotpe_formula1+num_isotope_formula2))*natural_abundance^(num_isotpe_formula1+num_isotope_formula2)

    ratio = ratio * (inten_formula2/inten_formula1)

  }

  return(ratio)
}

#' my_break_formula
#'
#' @param formula any chemical formula
#' @noRd
#' @return a list containing elem and count
#'
#' @examples my_break_formula("[13]C1C-1H4N2")
my_break_formula = function(formula = "[13]C1C-1H4N2"){
  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }

  return(list(elem=element2, count=number2))

}



#my function to process formula
{
  table_to_formula = function(temp_formula){
    r = c()
    elem = temp_formula$elem
    count = temp_formula$count
    len = length(elem)
    for(i in 1:len){
      if(count[i]!=0){
        r=paste(r, elem[i],  count[i], sep="")
      }
    }
    if(sum(count<0)!=0){
      r=paste(r, "Illegal_formula")
    }
    return (r)
  }
}


#Variables
{
  H_iso=(1.00782503224-1)
  N_iso=(14.0030740048-14)
  O_iso=(15.99491462-16)
  Cl_iso=(34.96885268-35)
  F_iso=(18.99840322-19)
  Br_iso=(78.9183371-79)
  I_iso=(126.904473-127)
  S_iso=(31.97207100-32)
  P_iso=(30.97376163-31)
  Si_iso=(27.97692653-28)
  B_iso=(11.009305167-11)
  Na_iso=(22.9897692809-23)
  K_iso=(38.96370668-39)
  Ca_iso=(39.962591-40)
  Ni_iso=(57.935343-58)
  Cu_iso=(62.929597-63)

  H_mass = 1.00782503224
  e_mass = 0.00054857990943


  # element ratio from 7 golden rule
  H_C_ratio_min = 0.1
  H_C_ratio_max = 6
  N_C_ratio_max = 4
  O_C_ratio_max = 3
  P_C_ratio_max = 2
  S_C_ratio_max = 3
  Si_C_ratio_max = 1
  F_C_ratio_max = 6
  Cl_C_ratio_max = 2
}


#' mz_formula
#'
#' @param Accurate_mass accurate mass from ms data
#' @param charge 1 for positive, 0 for neutral mass, -1 for negative
#' @param ppm output result tolerance window
#' @param C_range carbon number range
#' @param H_range hydrogen number range
#' @param O_range oxygen number range
#' @param N_range nitrogen number range
#' @param Cl_range chlorine number range
#' @param P_range phosphorus number range
#' @param S_range sulfur number range
#' @param Na_range sodium number range
#' @param K_range potassium number range
#' @param F_range fluorine number range
#' @param Br_range bromine number range
#' @param I_range iodine number range
#' @param Si_range silicon number range
#' @param N_rule nitrogen rule
#' @param db_min min of double bond and ring number
#' @param db_max max of double bond and ring number
#' @param B_range boron number range
#' @param metal_ion total number of metal ion
#' @param Ca_range Ca number range
#' @param Cu_range Cu number range
#' @param Ni_range Ni number range
#' @param Elem_ratio_rule element ratio from 7 golden rule, 0.1<= H/C <=6, N/C<=4, O/C<=3, P/C<=2, S/C<=3, Si/C<=1, F/C<=6, Cl/C<=2
#' @noRd
#' @return a dataframe recording formula, mass difference and db_r number that fall within ppm of input mass
#'
#' @importFrom dplyr bind_rows
#' @examples mz_formula(148.0604, 1, 5)

mz_formula = function(Accurate_mass = 148.0604,
                      charge = 1,
                      ppm = 5,
                      C_range = 0:99,
                      H_range = 0:100,
                      O_range = 0:20,
                      N_range = 0:12,
                      Cl_range = 0:0,
                      P_range = 0:3,
                      S_range = 0:3,
                      Na_range = 0:0,
                      K_range = 0:0,
                      F_range = 0:0,
                      Br_range = 0:0,
                      I_range = 0:0,
                      Si_range = 0:0,
                      B_range = 0:0,
                      Ca_range = 0:0,
                      Cu_range = 0:0,
                      Ni_range = 0:0,
                      N_rule = T,
                      Elem_ratio_rule = F,
                      db_min = 0,
                      db_max = 99,
                      metal_ion = 0:3)
{


  if(charge==0){mz_neutral=Accurate_mass}
  else{mz_neutral = (Accurate_mass - (H_mass-e_mass)*sign(charge))*abs(charge)}
  tolerance = mz_neutral*ppm/10^6

  H_min0 = min(H_range)
  C_min = min(C_range)
  C_max = max(C_range)

  temp_formula = list(elem = c("C", "H", "N", "O", "Cl", "S", "P",
                               "Na", "K", "F", "Br", "I", "Si","B",
                               "Ca","Cu","Ni"), count = rep(1,17))
  temp_formula$elem = temp_formula$elem[order(temp_formula$elem)]
  temp_formula_list = list()

  iteration=0
  mz_integer = floor(mz_neutral+0.2)
  mz_decimal = mz_neutral - mz_integer

  for(B in B_range){
    for(Si in Si_range){
      for(I in I_range){
        for(Br in Br_range){
          for(fluorine in F_range){
            for(Cu in Cu_range)
              for(Ni in Ni_range)
                for(Ca in Ca_range)
                  for(K in K_range){
                    for(Na in Na_range){
                      if((K + Na + Ca + Ni + Cu) > max(metal_ion)) next
                      CHNOPSCl = mz_integer - 23*Na - 39*K - 19*fluorine - 79*Br - 127*I - 28*Si - 11*B - 40*Ca - 58*Ni - 63*Cu
                      for(P in P_range){
                        for(S in S_range){
                          for(Cl in Cl_range){
                            CHNO = CHNOPSCl - 35*Cl - 32*S - 31*P
                            if(CHNO<0) break
                            for(N in N_range){
                              #(1) N rule
                              if(N_rule){
                                if(N%%2!=mz_integer%%2) next
                              }
                              for(O in O_range){
                                CH = CHNO - 14*N - 16*O

                                if(CH<0) break
                                else {
                                  Heavy_atom_decimal = (Cl_iso*Cl) + (F_iso*fluorine) + (Br_iso*Br) + (I_iso*I) + (S_iso*S) +
                                    (P_iso*P) + (Si_iso*Si) + (Na_iso*Na) + (K_iso*K) + (N_iso*N) + (O_iso*O) + (B_iso*B) +
                                    (Ca_iso*Ca) + (Cu_iso*Cu) + (Ni_iso*Ni)
                                  H_min = max(H_min0, floor((mz_decimal-Heavy_atom_decimal-tolerance)/H_iso))
                                  C_lower = max(C_min, floor(CH/14)-1)
                                  C_upper = min(C_max, floor(CH/12)+1)
                                }
                                for(C in C_lower:C_upper){

                                  iteration = iteration+1

                                  ##(2) H number
                                  H = CH - 12*C
                                  if(H < H_min) next

                                  ##(3) Saturation, ring&double bond
                                  db_r = (C+Si+1) - ((H+Cl+fluorine+Br+I+Na+K+Ca*2+Cu*2+Ni*2) - (B+N+P))/2;
                                  if(db_r < db_min | db_r >db_max) next

                                  ##(4) Decimal place
                                  differ = (H_iso*H) + Heavy_atom_decimal - mz_decimal
                                  if(differ > tolerance | differ < -tolerance) next

                                  ##(5) Element ratio from 7 golded rule
                                  if(Elem_ratio_rule & C!=0){
                                    Meet_ratio_rule = all (H/C >= H_C_ratio_min,
                                                           H/C <= H_C_ratio_max,
                                                           N/C <= N_C_ratio_max,
                                                           O/C <= O_C_ratio_max,
                                                           P/C <= P_C_ratio_max,
                                                           S/C <= S_C_ratio_max,
                                                           Si/C <= Si_C_ratio_max,
                                                           fluorine/C <= F_C_ratio_max,
                                                           Cl/C <=Cl_C_ratio_max
                                    )
                                    if(!Meet_ratio_rule) next
                                  }
                                  ## record

                                  temp_formula$count = c(B, Br, C, Ca, Cl, Cu, fluorine, H, I, K, N, Na, Ni, O, P, S, Si)
                                  temp_formula_list[[length(temp_formula_list)+1]]=c(table_to_formula(temp_formula),differ,db_r)
                                }}}}}}}}}}}}}

  if(mz_neutral > 700){
    mz_integer = floor(mz_neutral-.8)
    mz_decimal = mz_neutral - mz_integer
    for(B in B_range){
      for(Si in Si_range){
        for(I in I_range){
          for(Br in Br_range){
            for(fluorine in F_range){
              for(Cu in Cu_range)
                for(Ni in Ni_range)
                  for(Ca in Ca_range)
                    for(K in K_range){
                      for(Na in Na_range){
                        if((K + Na + Ca + Ni + Cu) > max(metal_ion)) next
                        CHNOPSCl = mz_integer - 23*Na - 39*K - 19*fluorine - 79*Br - 127*I - 28*Si - 11*B - 40*Ca - 58*Ni - 63*Cu
                        for(P in P_range){
                          for(S in S_range){
                            for(Cl in Cl_range){
                              CHNO = CHNOPSCl - 35*Cl - 32*S - 31*P
                              if(CHNO<0) break
                              for(N in N_range){
                                #(1) N rule
                                if(N_rule){
                                  if(N%%2!=mz_integer%%2) next
                                }
                                for(O in O_range){
                                  CH = CHNO - 14*N - 16*O

                                  if(CH<0) break
                                  else {
                                    Heavy_atom_decimal = (Cl_iso*Cl) + (F_iso*fluorine) + (Br_iso*Br) + (I_iso*I) + (S_iso*S) +
                                      (P_iso*P) + (Si_iso*Si) + (Na_iso*Na) + (K_iso*K) + (N_iso*N) + (O_iso*O) + (B_iso*B) +
                                      (Ca_iso*Ca) + (Cu_iso*Cu) + (Ni_iso*Ni)
                                    H_min = max(H_min0, floor((mz_decimal-Heavy_atom_decimal-tolerance)/H_iso))
                                    C_lower = max(C_min, floor(CH/14)-1)
                                    C_upper = min(C_max, floor(CH/12)+1)
                                  }
                                  for(C in C_lower:C_upper){

                                    iteration = iteration+1

                                    ##(2) H number
                                    H = CH - 12*C
                                    if(H < H_min) next

                                    ##(3) Saturation, ring&double bond
                                    db_r = (C+Si+1) - ((H+Cl+fluorine+Br+I+Na+K+Ca*2+Cu*2+Ni*2) - (B+N+P))/2;
                                    if(db_r < db_min | db_r >db_max) next

                                    ##(4) Decimal place
                                    differ = (H_iso*H) + Heavy_atom_decimal - mz_decimal
                                    if(differ > tolerance | differ < -tolerance) next

                                    ##(5) Element ratio from 7 golded rule
                                    if(Elem_ratio_rule & C!=0){
                                      Meet_ratio_rule = all (H/C >= H_C_ratio_min,
                                                             H/C <= H_C_ratio_max,
                                                             N/C <= N_C_ratio_max,
                                                             O/C <= O_C_ratio_max,
                                                             P/C <= P_C_ratio_max,
                                                             S/C <= S_C_ratio_max,
                                                             Si/C <= Si_C_ratio_max,
                                                             fluorine/C <= F_C_ratio_max,
                                                             Cl/C <=Cl_C_ratio_max
                                      )
                                      if(!Meet_ratio_rule) next
                                    }

                                    ## record

                                    temp_formula$count = c(B, Br, C, Ca, Cl, Cu, fluorine, H, I, K, N, Na, Ni, O, P, S, Si)
                                    temp_formula_list[[length(temp_formula_list)+1]]=c(table_to_formula(temp_formula),differ,db_r)
                                  }}}}}}}}}}}}}
  }





  #print(paste("iteration =", iteration))
  if(length(temp_formula_list)==0){
    return(data.frame(formula = as.character(), differ = as.numeric(), db_r = as.numeric()))
  }

  formula_df=as.data.frame(matrix(unlist(temp_formula_list),ncol=3,byrow = T),stringsAsFactors=F)
  colnames(formula_df) = c("formula", "differ", "db_r")
  formula_df$differ=as.numeric(formula_df$differ)
  formula_df$db_r=as.numeric(formula_df$db_r)
  formula_df = formula_df[with(formula_df, order(abs(differ))),]

  return(formula_df)
}

#' elem_num_query
#'
#' @param formula for example "C2H4O2"
#' @param elem_query it refers to element in query
#'
#' @return It returns the number of element in formula
#' @noRd
#'
#' @examples elem_num_query("C2H4O2","C")
elem_num_query = function(formula,elem_query){


  if(!is.character(formula)|is.na(formula)){return(NA)}

  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    #browser()
    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1
  }

  if(any(element2 == elem_query )){
    num_query = number2[elem_query == element2]

  }
  else{
    num_query=0

  }
  return(num_query)
}





