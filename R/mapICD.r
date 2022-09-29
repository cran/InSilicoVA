#' Map ICD-10 codes into the WHO 2016 cause list 
#' 
#' @param x a character object or a vector of ICD-10 codes  
#' @examples
#' mapICD("A90")
#' mapICD(c("A90", "C30"))
#' @export
#' 
mapICD <- function(x){

    icd10Map <- data.frame(
        rbind(
            c('Diarrhoeal diseases', 'A', 00, 09),
            c('Pulmonary tuberculosis', 'A', 15, 16),
            c('Other and unspecified infect dis', 'A', 17, 19),
            c('Other and unspecified infect dis', 'A', 20, 32),
            c('Tetanus (neonatorum)', 'A', 33, 33),
            c('Tetanus', 'A', 34, 35),
            c('Other and unspecified infect dis', 'A', 36, 36),
            c('Pertussis', 'A', 37, 37),
            c('Other and unspecified infect dis', 'A', 38, 38),
            c('Meningitis and encephalitis', 'A', 39, 39),
            c('Sepsis (non-obstetric)', 'A', 40, 41),
            c('Other and unspecified infect dis', 'A', 42, 89),
            c('Dengue fever', 'A', 90, 91),
            c('Haemorrhagic fever (non-dengue)', 'A', 92, 99),
            c('Other and unspecified infect dis', 'B', 00, 04),
            c('Measles', 'B', 05, 05),
            c('Other and unspecified infect dis', 'B', 06, 19),
            c('HIV/AIDS related death', 'B', 20, 24),
            c('Other and unspecified infect dis', 'B', 25, 49),
            c('Malaria', 'B', 50, 54),
            c('Other and unspecified infect dis', 'B', 55, 99),
            c('Oral neoplasms', 'C', 00, 06),
            c('Other and unspecified neoplasms', 'C', 07, 14),
            c('Digestive neoplasms', 'C', 15, 26),
            c('Respiratory neoplasms', 'C', 30, 39),
            c('Other and unspecified neoplasms', 'C', 40, 49),
            c('Breast neoplasms', 'C', 50, 50),
            c('Reproductive neoplasms (female)', 'C', 51, 58),
            c('Reproductive neoplasms (male)', 'C', 60, 63),
            c('Other and unspecified neoplasms', 'C', 60, 96),
            c('Other and unspecified neoplasms', 'D', 00, 48),
            c('Severe anaemia', 'D', 50, 56),
            c('Sickle cell with crisis', 'D', 57, 57),
            c('Severe anaemia', 'D', 58, 64),
            c('Other and unspecified NCD', 'D', 65, 89),
            c('Other and unspecified NCD', 'E', 00, 07),
            c('Diabetes mellitus', 'E', 10, 14),
            c('Other and unspecified NCD', 'E', 15, 35),
            c('Severe malnutrition', 'E', 40, 46),
            c('Other and unspecified NCD', 'E', 50, 90),
            c('Other and unspecified NCD', 'F', 00, 99),
            c('Meningitis and encephalitis', 'G', 00, 05),
            c('Other and unspecified NCD', 'G', 06, 37),
            c('Epilepsy', 'G', 40, 41),
            c('Other and unspecified NCD', 'G', 50, 99),
            c('Other and unspecified NCD', 'H', 00, 95),
            c('Other and unspecified cardiac dis', 'I', 00, 15),
            c('Acute cardiac disease', 'I', 20, 25),
            c('Other and unspecified cardiac dis', 'I', 26, 52),
            c('Stroke', 'I', 60, 69),
            c('Other and unspecified cardiac dis', 'I', 70, 99),
            c('Acute resp infect incl pneumonia', 'J', 00, 22),
            c('Other and unspecified NCD', 'J', 30, 39),
            c('Chronic obstructive pulmonary dis', 'J', 40, 44),
            c('Asthma', 'J', 45, 46),
            c('Other and unspecified NCD', 'J', 47, 99),
            c('Other and unspecified NCD', 'K', 00, 31),
            c('Other and unspecified NCD', 'K', 35, 38),
            c('Other and unspecified NCD', 'K', 40, 69),
            c('Liver cirrhosis', 'K', 70, 76),
            c('Other and unspecified NCD', 'K', 77, 93),
            c('Other and unspecified NCD', 'L', 00, 99),
            c('Other and unspecified NCD', 'M', 00, 99),
            c('Other and unspecified NCD', 'N', 00, 16),
            c('Renal failure', 'N', 17, 19),
            c('Other and unspecified NCD', 'N', 20, 99),
            c('Ectopic pregnancy', 'O', 00, 00),
            c('Other and unspecified maternal CoD', 'O', 01, 02),
            c('Abortion-related death', 'O', 03, 08),
            c('Pregnancy-induced hypertension', 'O', 10, 16),
            c('Other and unspecified maternal CoD', 'O', 20, 45),
            c('Obstetric haemorrhage', 'O', 46, 46),
            c('Other and unspecified maternal CoD', 'O', 47, 62),
            c('Obstructed labour', 'O', 63, 66),
            c('Obstetric haemorrhage', 'O', 67, 67),
            c('Other and unspecified maternal CoD', 'O', 68, 70),
            c('Ruptured uterus', 'O', 71, 71),
            c('Obstetric haemorrhage', 'O', 72, 72),
            c('Other and unspecified maternal CoD', 'O', 73, 75.2),
            c('Pregnancy-related sepsis', 'O', 75.3, 75.3),
            c('Other and unspecified maternal CoD', 'O', 75.4, 84),
            c('Pregnancy-related sepsis', 'O', 85, 85),
            c('Other and unspecified maternal CoD', 'O', 86, 98.9),
            c('Anaemia of pregnancy', 'O', 99.0, 99.0),
            c('Other and unspecified maternal CoD', 'O', 99.1, 99.8),
            c('Other and unspecified neonatal CoD', 'P', 00, 04),
            c('Prematurity', 'P', 05, 07),
            c('Other and unspecified neonatal CoD', 'P', 08, 15),
            c('Birth asphyxia', 'P', 20, 22),
            c('Neonatal pneumonia', 'P', 23, 25),
            c('Other and unspecified neonatal CoD', 'P', 26, 35),
            c('Neonatal sepsis', 'P', 36, 36),
            c('Other and unspecified neonatal CoD', 'P', 37, 94),
            c('Fresh stillbirth', 'P', 95, 95),
            c('Macerated stillbirth', 'P', 95, 95),  ## problem (same as above) CoD will always be Macerated (never Fresh)
            c('Other and unspecified neonatal CoD', 'P', 96, 96),
            c('Congenital malformation', 'Q', 00, 99),
            c('Other and unspecified NCD', 'R', 00, 09),
            c('Acute abdomen', 'R', 10, 10),
            c('Other and unspecified NCD', 'R', 11, 94),
            c('Unknown', 'R', 95, 99),
            c('Other and unspecified external CoD', 'S', 00, 99),
            c('Other and unspecified external CoD', 'T', 00, 99),
            c('Road traffic accident', 'V', 01, 89),
            c('Other transport accident', 'V', 90, 99),
            c('Accid fall', 'W', 00, 19),
            c('Other and unspecified external CoD', 'W', 20, 64),
            c('Accid drowning and submersion', 'W', 65, 74),
            c('Other and unspecified external CoD', 'W', 75, 99),
            c('Accid expos to smoke fire & flame', 'X', 00, 19),
            c('Contact with venomous plant/animal', 'X', 20, 29),
            c('Exposure to force of nature', 'X', 30, 39),
            c('Accid poisoning & noxious subs', 'X', 40, 49),
            c('Other and unspecified external CoD', 'X', 50, 59),
            c('Intentional self-harm', 'X', 60, 84),
            c('Assault', 'X', 85, 99),
            c('Assault', 'Y', 00, 09),
            c('Other and unspecified external CoD', 'Y', 10, 98)
        ), stringsAsFactors = FALSE
    )
    
    out <- rep(NA, length(x))
    for(jj in 1:length(x)){
        icd10Code <- x[jj]
        parsedCode <- cbind(substr(icd10Code, 1, 1),
                        substr(icd10Code, 2, 3),
                        substr(icd10Code, 5, 5))
        tmp <- 'Other (not in InterVA-5 cause list)'
        for (i in 1:nrow(icd10Map)) {  

            if (parsedCode[1] != icd10Map[i, 2]) next

            start <- as.numeric(icd10Map[i, 3])

            stop <- as.numeric(icd10Map[i, 4])
            if (round(stop) == stop) {
                stop <- as.numeric(paste0(stop, '.9'))
                }
            
            numCode <- paste0(parsedCode[2], '.', parsedCode[3])
            numCode <- as.numeric(numCode)
            if (start <= numCode & numCode <= stop) {
                if(icd10Map[i, 1] == "Tetanus (neonatorum)") icd10Map[i, 1] = "Tetanus"
                if(icd10Map[i, 1] == "Reproductive neoplasms (female)") icd10Map[i, 1] = "Reproductive neoplasms MF"
                if(icd10Map[i, 1] == "Reproductive neoplasms (male)") icd10Map[i, 1] = "Reproductive neoplasms MF"
                tmp <- icd10Map[i, 1]
                break
            }
        }
        out[jj] <- tmp
    }
    
    return(out)
}