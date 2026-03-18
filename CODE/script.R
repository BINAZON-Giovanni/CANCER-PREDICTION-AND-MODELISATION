#############################################################################################################################
#############################################################################################################################
                              
                                  # PROJET DE MODELISATION STATISTIQUE DE DONNEES COMPLEXES

#############################################################################################################################
#############################################################################################################################





######################################### Installation et importation des librairies  #######################################################

install.packages("labelled")

install.packages("rpart")

install.packages("rattle")

install.packages("pROc")

install.packages("adabag")

install.packages("ada")

install.packages("factoextra")

install.packages("missMDA")

library(ada)       
                          library(adabag)

library(pROC)        
                          library(rattle)

library(rpart)            
                          library(labelled)

library(dplyr)            
                          library(ggplot2)

library(tidyverse)
                          library(FactoMineR)

library(factoextra)
                          library(missMDA)

library(dplyr)
                          library(gt)

library(scales)

library(MASS)

############################################################################################################################
#                                 Prétraitement des données- description des données
############################################################################################################################



################################### Phase d'importation dataset et labellisation  #######################################################

cohorte = readRDS("Data/CohorteB")

colnames(cohorte)

str(cohorte)

cohorte$cible <- as.factor(cohorte$cible)


var_label(cohorte$RIBov) <- "Ratio de Pratique de Bovin"
var_label(cohorte$RIMou) <- "Ratio de Pratique de Moutons/Chèvre"
var_label(cohorte$RICoc) <- "Ratio de Pratique de Cochon"
var_label(cohorte$RIChe) <- "Ratio de Pratique de Chevaux"
var_label(cohorte$RIVol) <- "Ratio de Pratique de Volaille"
var_label(cohorte$RIPrai) <- "Ratio de Pratique de Prairies"
var_label(cohorte$RIVigne) <- "Ratio de Pratique de Vigne"
var_label(cohorte$RIMais) <- "Ratio de Pratique de Mais"
var_label(cohorte$RIBle) <- "Ratio de Pratique de Blé"
var_label(cohorte$RIPois) <- "Ratio de Pratique de Pois Fourragers/"
var_label(cohorte$RIBet) <- "Ratio de Pratique de Betteraves sucrières"
var_label(cohorte$RITou) <- "Ratio de Pratique de Tournesol"
var_label(cohorte$RICol) <- "Ratio de Pratique de Colza"
var_label(cohorte$RITabac) <- "Ratio de Pratique de Tabac"
var_label(cohorte$RIArb) <- "Ratio de Pratique de Arboriculture"
var_label(cohorte$RIPdt) <- "Ratio de Pratique de Pomme de Terre"
var_label(cohorte$RILegChamp) <- "Ratio de Pratique des autres cultures légumières en plein champs"
var_label(cohorte$RISerres) <- "Ratio de Pratique de la  cultures sous serres"
var_label(cohorte$cible) <- "Pésence ou Non de la maladie"



######################################### phase de nettoyage des données  #######################################################

summary(cohorte)

vars <- c("RIBov","RIMou","RICoc","RIChe","RIVol","RIPrai",
          "RIVigne","RIMais","RIBle","RIPois","RIBet","RITou",
          "RICol","RITabac","RIArb","RIPdt","RILegChamp","RISerres")


varUniQuanti <- function(variable, titreGraph = "Densité", titreX = "Variable",
                         show_hist = TRUE, show_mean = TRUE, show_median = TRUE,
                         col_fill = "#4C72B0", col_line = "#2A4E6E") {
  
    # Densité
    dens <- density(variable, na.rm = TRUE)
    
    # Histogramme optionnel
    if (show_hist) {
      hist(variable,
           freq = FALSE,
           col = adjustcolor(col_fill, alpha.f = 0.3),
           border = "white",
           main = titreGraph,
           xlab = titreX)
      lines(dens, lwd = 2, col = col_line)
    } else {
      plot(dens,
           main = titreGraph,
           xlab = titreX,
           lwd = 2,
           col = col_line)
    }
    
    # Ajout moyenne
    if (show_mean) {
      abline(v = mean(variable, na.rm = TRUE),
             col = "#D55E00", lwd = 2, lty = 2)
    }
    
    # Ajout médiane
    if (show_median) {
      abline(v = median(variable, na.rm = TRUE),
             col = "#009E73", lwd = 2, lty = 3)
    }
    
    legend("topright",
           legend = c("Moyenne", "Médiane"),
           col = c("#D55E00", "#009E73"),
           lty = c(2, 3),
           lwd = 2,
           bty = "n")
    
    # Quantiles
    q <- quantile(variable, probs = seq(0, 1, by = 0.01), na.rm = TRUE)
    return(q)
  }




export_outliers <- function(data, variable_name, seuil = 1,
                            dossier_sortie = "INFORMATION VALEURS ABBERANTES") {
  
  # Vérification de la variable
  if (!variable_name %in% names(data)) {
    stop(paste("La variable", variable_name, "n'existe pas dans le dataset."))
  }
  
  # Création du dossier si nécessaire
  if (!dir.exists(dossier_sortie)) {
    dir.create(dossier_sortie)
  }
  
  # Copie du dataset
  data_export <- data
  
  # Création de l'ID si absent
  if (!"ID" %in% names(data_export)) {
    data_export$ID <- rownames(data_export)
  }
  
  # Détection des valeurs aberrantes
  outliers <- data_export[data_export[[variable_name]] > seuil, ]
  
  # Sélection des colonnes voulues : ID + cible + variable analysée
  colonnes_a_garder <- c("ID", variable_name, "cible")
  outliers <- outliers[, colonnes_a_garder]
  
  # Nom du fichier
  nom_fichier <- paste0(dossier_sortie, "/Var_Abberants_", variable_name, ".csv")
  
  # Export CSV
  write.csv(outliers, nom_fichier, row.names = FALSE)
  
  # Création automatique d’un data.frame dans l’environnement global
  nom_objet <- paste0(variable_name, "_aberrants")
  assign(nom_objet, outliers, envir = .GlobalEnv)
  
  message("Objet créé dans R  ")
  
  
  message("Fichier exporté avec succès ")
  
  return(outliers)
}


# Nettoyage Ratio Bovin
varUniQuanti(cohorte$RIBov,
             "Courbe de densité et quantiles du ratio de la pratique de l'élevage de Bovin",
             "Variable Ratio Bovin")




# Nettoyage Ratio Mouton
varUniQuanti(cohorte$RIMou,
             "Courbe de densité et quantiles du ratio de la pratique de l'élevage de Mouton",
             "Variable Ratio Mouton")




# Nettoyage Ratio Cochon
varUniQuanti(cohorte$RICoc,
             "Courbe de densité et quantiles du ratio de la pratique de l'élevage de Cochon",
             "Variable Ratio Cochon")




# Nettoyage Ratio Cheval
varUniQuanti(cohorte$RIChe,
             "Courbe de densité et quantiles du ratio de la pratique de l'élevage de Chevaux",
             "Variable Ratio Chevaux")




# Nettoyage Ratio Volaille
varUniQuanti(cohorte$RIVol,
             "Courbe de densité et quantiles du ratio de la pratique de l'élevage de Volaille",
             "Variable Ratio Volaille")




# Nettoyage Ratio Prairie
varUniQuanti(cohorte$RIPrai,
             "Courbe de densité et quantiles du ratio de la pratique de la Prairie",
             "Variable Ratio Prairie")




# Nettoyage Ratio Vigne
varUniQuanti(cohorte$RIVigne,
             "Courbe de densité et quantiles du ratio de la pratique de la Vigne",
             "Variable Ratio Vigne")




# Nettoyage Ratio Maïs
varUniQuanti(cohorte$RIMais,
             "Courbe de densité et quantiles du ratio de la pratique du Maïs",
             "Variable Ratio Maïs")




# Nettoyage Ratio Blé
varUniQuanti(cohorte$RIBle,
             "Courbe de densité et quantiles du ratio de la pratique du Blé",
             "Variable Ratio Blé")




# Nettoyage Ratio Pois
varUniQuanti(cohorte$RIPois,
             "Courbe de densité et quantiles du ratio de la pratique du Pois",
             "Variable Ratio Pois")




# Nettoyage Ratio Betterave
varUniQuanti(cohorte$RIBet,
             "Courbe de densité et quantiles du ratio de la pratique de la Betterave",
             "Variable Ratio Betterave")




# Nettoyage Ratio Tournesol
varUniQuanti(cohorte$RITou,
             "Courbe de densité et quantiles du ratio de la pratique du Tournesol",
             "Variable Ratio Tournesol")




# Nettoyage Ratio Colza
varUniQuanti(cohorte$RICol,
             "Courbe de densité et quantiles du ratio de la pratique du Colza",
             "Variable Ratio Colza")




# Nettoyage Ratio Tabac
varUniQuanti(cohorte$RITabac,
             "Courbe de densité et quantiles du ratio de la pratique du Tabac",
             "Variable Ratio Tabac")




# Nettoyage Ratio Arboriculture
varUniQuanti(cohorte$RIArb,
             "Courbe de densité et quantiles du ratio de la pratique de l'Arboriculture",
             "Variable Ratio Arboriculture")




# Nettoyage Ratio Pomme de terre
varUniQuanti(cohorte$RIPdt,
             "Courbe de densité et quantiles du ratio de la pratique de la Pomme de Terre",
             "Variable Ratio Pomme de Terre")




# Nettoyage Ratio Légumes de plein champ
varUniQuanti(cohorte$RILegChamp,
             "Courbe de densité et quantiles du ratio de la pratique des Légumes de plein champ",
             "Variable Ratio Légumes plein champ")




# Nettoyage Ratio Serres
varUniQuanti(cohorte$RISerres,
             "Courbe de densité et quantiles du ratio de la pratique des Serres",
             "Variable Ratio Serres")






vars_acp <- cohorte[, c("RIBov","RIMou","RICoc","RIChe","RIVol","RIPrai",
                        "RIVigne","RIMais","RIBle","RIPois","RIBet","RITou",
                        "RICol","RITabac","RIArb","RIPdt","RILegChamp","RISerres")]



res_acp <- PCA(vars_acp, 
               scale.unit = TRUE,   # Standardisation
               ncp = 5,             # Nombre de composantes à calculer
               graph = FALSE)



fviz_pca_var(res_acp,
             col.var = "cos2", # Coloration selon qualité de représentation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)






lapply(vars, function(v) export_outliers(cohorte, v))


######################################### CREATION DU TABLEAU DE SYNTHESE  #######################################################

 
# Variables quantitatives
vars_quanti <- c("RIBov","RIMou","RICoc","RIChe","RIVol","RIPrai",
                 "RIVigne","RIMais","RIBle","RIPois","RIBet","RITou",
                 "RICol","RITabac","RIArb","RIPdt","RILegChamp","RISerres")

# Fonction de résumé univarié
resume_univarie <- function(data, vars) {
  
  resultat <- lapply(vars, function(v) {
    
    x <- data[[v]]
    
    Effectif <- sum(!is.na(x))
    Nb_aberrantes <- sum(x > 1, na.rm = TRUE)
    Pourcentage_aberrantes <- (Nb_aberrantes / Effectif) * 100
    Pourcentage_non_nuls <- mean(x > 0, na.rm = TRUE) * 100
    
    data.frame(
      Variable = v,
      Effectif = Effectif,
      Moyenne = mean(x, na.rm = TRUE),
      Ecart_type = sd(x, na.rm = TRUE),
      Mediane = median(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE),
      Non_nuls = Pourcentage_non_nuls,
      Nb_aberrantes = Nb_aberrantes,
      Pourcentage_aberrantes = Pourcentage_aberrantes
    )
  })
  
  bind_rows(resultat)
}




# Création du tableau
table_univarie <- resume_univarie(cohorte, vars_quanti)

labels_variables <- c(
  RIBov = "Bovin",
  RIMou = "Moutons / Chèvres",
  RICoc = "Cochon",
  RIChe = "Chevaux",
  RIVol = "Volaille",
  RIPrai = "Prairies",
  RIVigne = "Vigne",
  RIMais = "Maïs",
  RIBle = "Blé",
  RIPois = "Pois fourragers",
  RIBet = "Betteraves sucrières",
  RITou = "Tournesol",
  RICol = "Colza",
  RITabac = "Tabac",
  RIArb = "Arboriculture",
  RIPdt = "Pomme de terre",
  RILegChamp = "Légumes plein champ",
  RISerres = "Cultures sous serres"
)

table_univarie$Variable <- labels_variables[table_univarie$Variable]


table_univarie %>%
  mutate(
    `Moyenne (Ecart_type)` = sprintf("%.4f (%.4f)", Moyenne, Ecart_type),
    Non_nuls = round(Non_nuls, 1),
    Pourcentage_aberrantes = round(Pourcentage_aberrantes, 2)
  ) %>%
  dplyr::select(
    Variable,
    Effectif,
    `Moyenne (Ecart_type)`,
    Mediane,
    Min,
    Max,
    Non_nuls,
    Nb_aberrantes,
    Pourcentage_aberrantes
  ) %>%
  
  gt() %>%
  tab_header(
    title = md("**Description univariée des variables quantitatives**"),
    subtitle = md("Ratios de pratiques agricoles")
  ) %>%
  
  cols_label(
    Non_nuls = "% non nuls",
    Pourcentage_aberrantes = "% aberrantes"
  ) %>%
  
  fmt_number(
    columns = c(Mediane, Min, Max),
    decimals = 4
  ) %>%
  
  data_color(
    columns = c(Non_nuls, Pourcentage_aberrantes),
    colors = col_numeric(
      palette = c("#E8F5E9", "#FFCDD2"),
      domain = NULL
    )
  ) %>%
  
  tab_source_note(
    source_note = md(
      "**Lecture du tableau :**  
      *% non nuls* correspond à la proportion d’observations strictement supérieures à 0.  
      *% aberrantes* correspond à la proportion de ratios strictement supérieurs à 1."
    )
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  
  opt_table_outline() %>%
  opt_all_caps()


################################## TRANSFORMATION DES VALEURS ABBERANTES EN MANQUANTES ####################################################


# Fonction de Modification des valeurs abbérantes par des NA puis 

replace_outliers_with_NA <- function(data, variable_name, seuil = 1) {
  
  
  if (!variable_name %in% names(data)) {
    stop(paste("La variable", variable_name, "n'existe pas dans le dataset."))
  }
  
  
  data_modif <- data
  
 
  data_modif[[variable_name]][data_modif[[variable_name]] > seuil] <- NA
  
  message("Valeurs aberrantes remplacées par NA pour : ", variable_name)
  
  return(data_modif)
}

for (v in vars) {
  cohorte <- replace_outliers_with_NA(cohorte, v)
}

# Fonction d'imputation des valeurs manquantes via l' ACP

impute_by_acp <- function(data, vars_quanti, ncp = 5) {
  
  # Vérification
  if (!all(vars_quanti %in% names(data))) {
    stop("Certaines variables quantitatives n'existent pas dans le dataset.")
  }
  
  # Extraction des variables quantitatives
  data_quanti <- data[, vars_quanti]
  
  # Imputation via ACP
  res_impute <- imputePCA(data_quanti, ncp = ncp)
  
  # Remplacement dans le dataset original
  data_imputed <- data
  data_imputed[, vars_quanti] <- res_impute$completeObs
  
  message("Imputation ACP réalisée avec succès.")
  
  return(data_imputed)
}

cohorte <- impute_by_acp(cohorte, vars, ncp = 5)



#####################################################################################################################################
######################################### PHASE DE PRE-SELECTION DE VARIABLE  #######################################################
#####################################################################################################################################



bootstrap_variable_selection <- function(data, cible, variables, 
                                         n_boot = 50, seuil_selection = 0.7) {
  
  # Stockage des sélections
  selections <- matrix(0, nrow = n_boot, ncol = length(variables))
  colnames(selections) <- variables
  
  for (i in 1:n_boot) {
    
    # Bootstrap : échantillon avec remise
    indices <- sample(1:nrow(data), replace = TRUE)
    data_boot <- data[indices, ]
    
    # Formule complète
    formule <- as.formula(
      paste(cible, "~", paste(variables, collapse = " + "))
    )
    
    # Modèle complet
    mod_full <- glm(formule, data = data_boot, family = binomial)
    
    # Sélection automatique stepwise
    mod_step <- stepAIC(mod_full, direction = "both", trace = FALSE)
    
    # Variables sélectionnées dans ce bootstrap
    vars_sel <- names(coef(mod_step))[-1]  # enlever l'intercept
    
    # Marquer les variables sélectionnées
    selections[i, vars_sel] <- 1
  }
  
  # Calcul des fréquences de sélection
  freq <- colMeans(selections)
  
  # Variables robustes
  vars_robustes <- names(freq[freq >= seuil_selection])
  
  # Résultat
  list(
    frequences = freq,
    variables_robustes = vars_robustes,
    tableau_complet = selections
  )
}


resultats_selection <- bootstrap_variable_selection(
  data = cohorte,
  cible = "cible",
  variables = vars,
  n_boot = 50,
  seuil_selection = 0.7
)

FrequencesVariables=data.frame(resultats_selection$frequences) 

vars_robustes <- resultats_selection$variables_robustes

vars <-vars_robustes

colonnes <- c("cible", resultats_selection$variables_robustes)

# Extraction du sous-dataset
dataset <- cohorte[, colonnes]


#################################################################################################################
############################## PHASE DE MODELISATION PAR REGRESSION LOGISTIQUE #############################################
#################################################################################################################




create_splits <- function(data, cible, n = 40, prop_train = 0.7) {
  
  splits <- list()
  
  for (i in 1:n) {
    set.seed(120+i)
    train_index <- sample(1:nrow(data), size = floor(prop_train * nrow(data)))
    
    splits[[i]] <- list(
      train = data[train_index, ],
      test  = data[-train_index, ]
    )
  }
  
  return(splits)
}



compute_metrics <- function(pred, truth) {
  
  pred  <- as.character(pred)
  truth <- as.character(truth)
  
  tab <- table(pred, truth)
  
  sens <- tab["1","1"] / sum(tab[,"1"])
  spec <- tab["0","0"] / sum(tab[,"0"])
  
  return(list(sens = sens, spec = spec))
}

############################# FUNCTION DE MODELISATION PAR REGRESSION LOGISTIQUE CLASSIQUE ######################################

run_logistic <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- glm(formule, data = train, family = binomial)
  
  prob <- predict(mod, newdata = test, type = "response")
  pred <- ifelse(prob > 0.5, 1, 0)
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}

############################# FONTION DE MODELISATION PAR REGRESSION LOGISTIQUE BAGGING ######################################
library(adabag)

run_bagging <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- bagging(formule, data = train)
  
  pred <- predict(mod, newdata = test)$class
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}

############################# FONCTION DE MODELISATION PAR REGRESSION LOGISTIQUE BOOSTING ######################################

run_boosting <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- boosting(formule, data = train)
  
  pred <- predict(mod, newdata = test)$class
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}

############################# FONCTION DE MODELISATION PAR REGRESSION LOGISTIQUE PRE - CLASSIFICATION ######################################
install.packages("NbClust")
install.packages("FNN")
library(FNN)
library(NbClust)

run_clustering <- function(train, vars, k = 3) {
  
  scaled <- scale(train[, vars])
  km <- kmeans(scaled, centers = k, nstart = 20)
  
  train$cluster <- km$cluster
  list(train = train, model = km)
}

run_logistic_by_cluster <- function(train, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  split(train, train$cluster) |>
    lapply(function(df) glm(formule, data = df, family = binomial))
}

choose_k_2_7 <- function(train, vars) {
  
  set.seed(123)
  idx <- sample(seq_len(nrow(train)), size = 0.2 * nrow(train))
  
  train_main <- train[-idx, ]
  valid      <- train[idx, ]
  
  train_scaled <- scale(train_main[, vars])
  valid_scaled <- scale(valid[, vars],
                        center = attr(train_scaled, "scaled:center"),
                        scale  = attr(train_scaled, "scaled:scale"))
  
  ks <- 2:7
  errors <- sapply(ks, function(k) {
    pred <- FNN::knn(train_scaled, valid_scaled, cl = train_main$cluster, k = k)
    mean(pred != valid$cluster)
  })
  
  ks[which.min(errors)]
}

assign_clusters_knn <- function(train, test, vars) {
  
  k <- choose_k_2_7(train, vars)
  message("k optimal choisi = ", k)
  
  train_scaled <- scale(train[, vars])
  test_scaled  <- scale(test[, vars],
                        center = attr(train_scaled, "scaled:center"),
                        scale  = attr(train_scaled, "scaled:scale"))
  
  pred <- tryCatch(
    FNN::knn(train_scaled, test_scaled, cl = train$cluster, k = k),
    error = function(e) {
      message("⚠️ KNN a échoué → cluster 1 assigné")
      rep(1, nrow(test))
    }
  )
  
  test$cluster <- as.numeric(pred)
  test
}

predict_by_cluster <- function(test, models) {
  
  sapply(seq_len(nrow(test)), function(i) {
    cl <- as.character(test$cluster[i])
    if (!cl %in% names(models)) cl <- names(models)[1]
    
    prob <- predict(models[[cl]], newdata = test[i, ], type = "response")
    as.numeric(prob > 0.5)
  })
}



run_preclassification <- function(train, test, cible, vars, k_clusters = 3) {
  
  # 1. Clustering du train
  clust <- run_clustering(train, vars, k = k_clusters)
  train <- clust$train
  
  # 2. Modèles logistiques par cluster
  models <- run_logistic_by_cluster(train, cible, vars)
  
  # 3. Affectation des clusters au test
  test <- assign_clusters_knn(train, test, vars)
  
  # 4. Prédictions
  preds <- predict_by_cluster(test, models)
  
  # 5. Métriques
  metrics <- compute_metrics(preds, test[[cible]])
  
  return(metrics)
}

############################## FONCTION PRINCIPALE DE MODELiSATION ##############################################################

run_all_models <- function(data, cible, vars, n = 40) {
  
  splits <- create_splits(data, cible, n)
  
  results <- data.frame(
    iteration = 1:n,
    sens_log = NA, spec_log = NA,
    sens_bag = NA, spec_bag = NA,
    sens_boost = NA, spec_boost = NA,
    sens_preclass = NA, spec_preclass = NA
  )
  
  for (i in 1:n) {
    
    train <- splits[[i]]$train
    test  <- splits[[i]]$test
    
    # Logistique
    log_res <- run_logistic(train, test, cible, vars)
    results$sens_log[i] <- log_res$sens
    results$spec_log[i] <- log_res$spec
    
    # Bagging
    bag_res <- run_bagging(train, test, cible, vars)
    results$sens_bag[i] <- bag_res$sens
    results$spec_bag[i] <- bag_res$spec
    
    # Boosting
    boost_res <- run_boosting(train, test, cible, vars)
    results$sens_boost[i] <- boost_res$sens
    results$spec_boost[i] <- boost_res$spec
    
    # Pré-classification
    pre_res <- run_preclassification(train, test, cible, vars)
    results$sens_preclass[i] <- pre_res$sens
    results$spec_preclass[i] <- pre_res$spec
    
  }
  
  return(results)
}


dataset = cohorte[1:10000, ]

resultats <- run_all_models(cohorte, "cible", vars_robustes, n = 40)


library(tidyverse)

# Transformation en format long
resultats_long <- resultats %>%
  pivot_longer(
    cols = -iteration,
    names_to = c("metrique", "modele"),
    names_sep = "_",
    values_to = "valeur"
  )

ggplot(resultats_long %>% filter(metrique == "sens"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparaison des sensibilités",
       x = "Méthode",
       y = "Sensibilité") +
  theme_minimal()

ggplot(resultats_long %>% filter(metrique == "spec"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparaison des spécificités",
       x = "Méthode",
       y = "Spécificité") +
  theme_minimal()




















#################################################################################################################
############################## PHASE DE MODELISATION PAR ARBRE DE DECISION #############################################
#################################################################################################################

############################# FUNCTION DE MODELISATION PAR ARBRE DE DECISION CLASSIQUE ######################################

library(rpart)

run_tree <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- rpart(formule, data = train, method = "class")
  
  pred <- predict(mod, newdata = test, type = "class")
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}


############################# FUNCTION DE MODELISATION PAR ARBRE DE DECISION BOOSTING ######################################


library(adabag)

run_boosting_tree <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- boosting(formule, data = train, mfinal = 50)
  
  pred <- predict(mod, newdata = test)$class
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}

############################# FUNCTION DE MODELISATION PAR ARBRE DE DECISION: FORET ALEATOIRE ######################################

install.packages("randomForest")
library(randomForest)

run_random_forest <- function(train, test, cible, vars) {
  
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  mod <- randomForest(formule, data = train, ntree = 300)
  
  pred <- predict(mod, newdata = test)
  
  metrics <- compute_metrics(pred, test[[cible]])
  
  return(metrics)
}

############################## FONCTION PRINCIPALE DE MODELiSATION ##############################################################

run_all_trees <- function(data, cible, vars, n = 40) {
  
  splits <- create_splits(data, cible, n)
  
  results <- data.frame(
    iteration = 1:n,
    sens_tree = NA, spec_tree = NA,
    sens_boost = NA, spec_boost = NA,
    sens_rf = NA, spec_rf = NA
  )
  
  for (i in 1:n) {
    
    train <- splits[[i]]$train
    test  <- splits[[i]]$test
    
    # Arbre classique
    tree_res <- run_tree(train, test, cible, vars)
    results$sens_tree[i] <- tree_res$sens
    results$spec_tree[i] <- tree_res$spec
    
    # Boosting
    boost_res <- run_boosting_tree(train, test, cible, vars)
    results$sens_boost[i] <- boost_res$sens
    results$spec_boost[i] <- boost_res$spec
    
    # Forêt aléatoire
    rf_res <- run_random_forest(train, test, cible, vars)
    results$sens_rf[i] <- rf_res$sens
    results$spec_rf[i] <- rf_res$spec
  }
  
  return(results)
}

resultats_trees <- run_all_trees(cohorte, "cible", vars_robustes, n = 40)

library(tidyverse)

resultats_trees_long <- resultats_trees %>%
  pivot_longer(
    cols = -iteration,
    names_to = c("metrique", "modele"),
    names_sep = "_",
    values_to = "valeur"
  )

ggplot(resultats_trees_long %>% filter(metrique == "sens"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparaison des sensibilités (arbres)",
       x = "Méthode",
       y = "Sensibilité") +
  theme_minimal()

ggplot(resultats_trees_long %>% filter(metrique == "spec"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparaison des spécificités (arbres)",
       x = "Méthode",
       y = "Spécificité") +
  theme_minimal()




#################################################################################################################
############################## PHASE DE MODELISATION PAR RESEAU DE NEURONNE #############################################
#################################################################################################################



scale_data <- function(train, test, vars) {
  
  means <- apply(train[, vars], 2, mean)
  sds   <- apply(train[, vars], 2, sd)
  
  train_scaled <- train
  test_scaled  <- test
  
  train_scaled[, vars] <- scale(train[, vars], center = means, scale = sds)
  test_scaled[, vars]  <- scale(test[, vars], center = means, scale = sds)
  
  return(list(train = train_scaled, test = test_scaled))
}

library(neuralnet)


architectures <- list(
  c(5),        # 1 couche cachée, 5 neurones
  c(10),       # 1 couche cachée, 10 neurones
  c(3, 3),     # 2 couches cachées : 5 puis 3 neurones
  c(3, 2)     # 2 couches cachées : 10 puis 5 neurones
)

run_neuralnet <- function(train, test, cible, vars, hidden_layers) {
  
  # Cible en 0/1 numérique
  train[[cible]] <- ifelse(train[[cible]] == 1, 1, 0)
  test[[cible]]  <- ifelse(test[[cible]] == 1, 1, 0)
  
  # Mise à l’échelle
  scaled <- scale_data(train, test, vars)
  train <- scaled$train
  test  <- scaled$test
  
  # Formule
  formule <- as.formula(paste(cible, "~", paste(vars, collapse = "+")))
  
  # Entraînement du réseau
  nn <- neuralnet(
    formule,
    data = train,
    hidden = hidden_layers,
    linear.output = FALSE,
    lifesign = "minimal",
    stepmax = 1e7,       # <-- IMPORTANT
    threshold = 0.05     # <-- IMPORTANT
    
  )
  
  # 🚨 Si le réseau n’a pas convergé → pas de prédiction possible
  if (is.null(nn$weights)) {
    warning("Le réseau n'a pas convergé : architecture trop complexe.")
    return(list(sens = NA, spec = NA))
  }
  
  
  # Prédiction
  prob <- compute(nn, test[, vars])$net.result
  
  # Si 2 colonnes → on garde la colonne 1
  if (is.matrix(prob) && ncol(prob) == 2) {
    prob <- prob[, 1]
  }
  
  prob <- as.vector(prob)
  pred <- ifelse(prob > 0.5, 1, 0)
  
  # Vérification
  if (length(pred) != length(test[[cible]])) {
    stop("Erreur : pred et truth n'ont pas la même longueur.")
  }
  
  metrics <- compute_metrics(pred, test[[cible]])
  return(metrics)
}


run_all_neuralnets <- function(data, cible, vars, architectures, n = 40) {
  
  splits <- create_splits(data, cible, n)
  
  results <- data.frame(
    iteration = 1:n
  )
  
  # Colonnes dynamiques selon le nombre d’architectures
  for (a in 1:length(architectures)) {
    results[[paste0("sens_nn", a)]] <- NA
    results[[paste0("spec_nn", a)]] <- NA
  }
  
  for (i in 1:n) {
    
    train <- splits[[i]]$train
    test  <- splits[[i]]$test
    
    # Pour chaque architecture
    for (a in 1:length(architectures)) {
      
      nn_res <- run_neuralnet(train, test, cible, vars, architectures[[a]])
      
      results[i, paste0("sens_nn", a)] <- nn_res$sens
      results[i, paste0("spec_nn", a)] <- nn_res$spec
    }
  }
  
  return(results)
}


resultats_nn <- run_all_neuralnets(
  data = cohorte,
  cible = "cible",
  vars = vars_robustes,
  architectures = architectures,
  n = 40
)


resultats_nn_long <- resultats_nn %>%
  pivot_longer(
    cols = -iteration,
    names_to = c("metrique", "modele"),
    names_sep = "_",
    values_to = "valeur"
  )

ggplot(resultats_nn_long %>% filter(metrique == "sens"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Sensibilité des réseaux de neurones",
       x = "Architecture",
       y = "Sensibilité") +
  theme_minimal()

ggplot(resultats_nn_long %>% filter(metrique == "spec"),
       aes(x = modele, y = valeur, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Spécificité des réseaux de neurones",
       x = "Architecture",
       y = "Spécificité") +
  theme_minimal() 
