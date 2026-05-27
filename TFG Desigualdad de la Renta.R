################################################################################
# TFG: SCRIPT Desigualdad de la renta equivalente en España (2016–2021):
# medición y evolución con microdatos del Panel de Hogares del IEF
################################################################################

library(tidyverse)
library(survey)
library(convey)
library(ineq)
library(pROC)

# --- 1. CONFIGURACIÓN ---
ruta_base <- "RUTA_A_LOS_DATOS"
# Modificar la variable "ruta_base" indicando la ubicación local
# donde se encuentran los microdatos del IEF/MCVL.

anios_estudio <- c(2016, 2017, 2018, 2019, 2020, 2021)

# ################################################################################
# --- BLOQUE 1: FUNCIÓN DE CARGA ---
# ################################################################################

procesar_limpio <- function(val_anio) {
  tryCatch({
    folder <- file.path(ruta_base, val_anio)
    f_iden <- list.files(folder, pattern = "IDEN", full.names = TRUE)[1]
    f_rent <- list.files(folder, pattern = "Renta", full.names = TRUE)[1]
    
    anchos_iden <- if(val_anio <= 2020) {
      fwf_cols(IDENPER=c(1,11), IDENHOG=c(12,35), FACTOR=c(36,55), ANNAC=c(76,79))
    } else {
      fwf_cols(IDENPER=c(1,11), IDENHOG=c(12,35), FACTOR=c(56,75), ANNAC=c(116,119))
    }
    
    p6 <- if(val_anio <= 2017) 467 else if(val_anio <= 2019) 503 else 539
    p7 <- p6 + 12
    
    anchos_rent <- fwf_cols(IDENPER=c(1,11), M1=c(119,130), M2=c(131,142), M3=c(191,202),
                            M4=c(227,238), M5=c(287,298), M6=c(p6,p6+11), M7=c(p7,p7+11))
    
    df_iden <- read_fwf(f_iden, anchos_iden, col_types = cols(IDENPER = col_character(), IDENHOG = col_character()), show_col_types=F)
    df_rent <- read_fwf(f_rent, anchos_rent, col_types = cols(IDENPER = col_character()), na=c(""," ","."), show_col_types=F)
    
    df_res <- df_iden %>%
      inner_join(df_rent, by = "IDENPER") %>%
      mutate(across(starts_with("M"), ~coalesce(as.numeric(.x), 0)),
             edad = val_anio - as.numeric(str_extract(ANNAC, "\\d+"))) %>%
      group_by(IDENHOG) %>%
      summarise(
        IDENPER = first(IDENPER),
        renta_hogar = sum((M1 + M2 + M3 + M4 + M5 + M6) - M7, na.rm=T),
        adultos_hogar = sum(edad >= 14, na.rm=T),
        ninos_hogar   = sum(edad < 14, na.rm=T),
        f_elev = as.numeric(first(FACTOR)), 
        .groups = "drop"
      ) %>%
      mutate(peso_ocde = 1 + 0.5 * pmax(0, adultos_hogar - 1) + 0.3 * ninos_hogar,
             renta_equiv = renta_hogar / peso_ocde, anio = val_anio) %>%
      filter(renta_equiv > 0 & !is.na(f_elev))
    
    return(df_res)
  }, error = function(e) return(NULL))
}


# Carga
lista_hogares <- map(anios_estudio, ~procesar_limpio(.x))
base_total <- bind_rows(lista_hogares)

# --- AUDITORÍA DE PERSISTENCIA ---
cat("\n--- AUDITORÍA DE SEGUIMIENTO (HOGARES) ---\n")
auditoria <- map_df(2:length(lista_hogares), function(i) {
  tibble(Periodo = paste(unique(lista_hogares[[i-1]]$anio), "-", unique(lista_hogares[[i]]$anio)),
         Hogares_Comunes = length(intersect(lista_hogares[[i-1]]$IDENHOG, lista_hogares[[i]]$IDENHOG)))
})
print(auditoria)

#  SANEAMIENTO (winsorización) ---
base_analisis <- base_total %>%
  group_by(anio) %>%
  filter(renta_equiv > quantile(renta_equiv, 0.01) & renta_equiv < quantile(renta_equiv, 0.99)) %>%
  ungroup()


# ################################################################################
# --- ÍNDICES DE GINI Y THEIL,INDICADORES PALMA RATIO Y S80/S20 
#     Y DESCOMPOSICIÓN HISTÓRICA CCAA) ---
# ################################################################################

# Diseño
diseno_evol <- convey_prep(svydesign(ids=~1, weights=~f_elev, data=base_analisis))

# ==============================================================================
# A. CÁLCULO DE INDICADORES GLOBALES DE LA SERIE HISTÓRICA (2016-2021)
# ==============================================================================
cat("\n--- EVOLUCIÓN HISTÓRICA DEL ÍNDICE DE GINI ---\n")
print(svyby(~renta_equiv, by=~anio, design=diseno_evol, FUN=svygini))

cat("\n--- EVOLUCIÓN HISTÓRICA DEL ÍNDICE DE THEIL (GE 1) ---\n")
print(svyby(~renta_equiv, by=~anio, design=diseno_evol, FUN=svygei, epsilon = 1))

cat("\n--- EVOLUCIÓN HISTÓRICA DEL ÍNDICE DE PALMA (P90/P40) ---\n")
calcular_palma_nativo <- function(df_sub) {
  p40 <- Hmisc::wtd.quantile(df_sub$renta_equiv, weights = df_sub$f_elev, probs = 0.40)
  p90 <- Hmisc::wtd.quantile(df_sub$renta_equiv, weights = df_sub$f_elev, probs = 0.90)
  
  renta_p90_rico <- sum(df_sub$renta_equiv[df_sub$renta_equiv >= p90] * df_sub$f_elev[df_sub$renta_equiv >= p90], na.rm = TRUE)
  renta_p40_pobre <- sum(df_sub$renta_equiv[df_sub$renta_equiv <= p40] * df_sub$f_elev[df_sub$renta_equiv <= p40], na.rm = TRUE)
  
  return(round(renta_p90_rico / renta_p40_pobre, 4))
}

tabla_palma_historica <- base_analisis %>%
  group_by(anio) %>%
  summarise(Palma_Ratio = calcular_palma_nativo(cur_data()), .groups = "drop")

print(as.data.frame(tabla_palma_historica))

cat("\n--- EVOLUCIÓN HISTÓRICA DEL ÍNDICE S80/S20 ---\n")
calcular_s80s20_nativo <- function(df_sub) {
  p20 <- Hmisc::wtd.quantile(df_sub$renta_equiv,weights = df_sub$f_elev,probs = 0.20,na.rm = TRUE)
  p80 <- Hmisc::wtd.quantile(df_sub$renta_equiv,weights = df_sub$f_elev,probs = 0.80, na.rm = TRUE  )
  renta_p20 <- sum( df_sub$renta_equiv[df_sub$renta_equiv <= p20] * df_sub$f_elev[df_sub$renta_equiv <= p20], na.rm = TRUE )
  renta_p80 <- sum(    df_sub$renta_equiv[df_sub$renta_equiv >= p80] *df_sub$f_elev[df_sub$renta_equiv >= p80],na.rm = TRUE)
  
  return(round(renta_p80 / renta_p20, 4))
}

tabla_s80s20 <- base_analisis %>%
  group_by(anio) %>%
  summarise( S80_S20 = calcular_s80s20_nativo(cur_data()), .groups = "drop" )

print(as.data.frame(tabla_s80s20))

# ==============================================================================
# B. DESCOMPOSICIÓN DINÁMICA DE THEIL POR CCAA, SERIE (2016-2021)
# ==============================================================================
cat("\n==================================================================\n")
cat("    DESCOMPOSICIÓN HISTÓRICA DE THEIL POR CCAA (SERIE 2016-2021)  \n")
cat("==================================================================\n")

# Función para calcular el Índice de Theil Ponderado 
calcular_theil_nativo <- function(x, w) {
  ind <- !is.na(x) & !is.na(w) & x > 0 & w > 0
  x <- x[ind]
  w <- w[ind]
  
  w_prop <- w / sum(w)
  media_w <- sum(x * w_prop)
  r_relativa <- x / media_w
  
  # Fórmula oficial de Entropía GE(1) / Theil
  theil <- sum(w_prop * r_relativa * log(r_relativa))
  return(theil)
}

ejecutar_descomposicion_anual <- function(val_anio) {
  file_raw_anio <- list.files(file.path(ruta_base, val_anio), pattern = "IDEN", full.names = TRUE)[1]
  
  df_ccaa_anio <- read_lines(file_raw_anio) %>%
    tibble(linea = .) %>%
    mutate(
      IDENPER = as.character(as.integer(str_trim(substr(linea, 1, 11)))),
      provincia_cod = substr(linea, 23, 24),
      CCAA = case_when(
        provincia_cod %in% c("04", "11", "14", "18", "21", "23", "29", "41") ~ "Andalucía",
        provincia_cod %in% c("22", "44", "50")                             ~ "Aragón",
        provincia_cod %in% c("33")                                         ~ "Asturias",
        provincia_cod %in% c("07")                                         ~ "Baleares",
        provincia_cod %in% c("35", "38")                                   ~ "Canarias",
        provincia_cod %in% c("39")                                         ~ "Cantabria",
        provincia_cod %in% c("05", "09", "24", "34", "37", "40", "42", "47", "49") ~ "Castilla y León",
        provincia_cod %in% c("02", "13", "16", "19", "45")                 ~ "Castilla-La Mancha",
        provincia_cod %in% c("08", "17", "25", "43")                       ~ "Cataluña",
        provincia_cod %in% c("03", "12", "46")                             ~ "C. Valenciana",
        provincia_cod %in% c("06", "10")                                   ~ "Extremadura",
        provincia_cod %in% c("15", "27", "32", "36")                       ~ "Galicia",
        provincia_cod %in% c("28")                                         ~ "Madrid",
        provincia_cod %in% c("30")                                         ~ "Murcia",
        provincia_cod %in% c("31")                                         ~ "Navarra",
        provincia_cod %in% c("01", "20", "48")                             ~ "País Vasco",
        provincia_cod %in% c("26")                                         ~ "La Rioja",
        TRUE                                                               ~ "Otras"
      )
    ) %>%
    distinct(IDENPER, .keep_all = TRUE) %>%
    select(IDENPER, CCAA)
  
  base_descomp_anio <- base_analisis %>%
    filter(anio == val_anio) %>%
    mutate(IDENPER = as.character(as.integer(str_trim(IDENPER)))) %>%
    inner_join(df_ccaa_anio, by = "IDENPER")
  
  if(nrow(base_descomp_anio) > 0) {
    theil_total <- calcular_theil_nativo(base_descomp_anio$renta_equiv, base_descomp_anio$f_elev)
    
    tabla_regiones <- base_descomp_anio %>%
      group_by(CCAA) %>%
      summarise(
        renta_tot = sum(renta_equiv * f_elev, na.rm = TRUE),
        Theil_Interno = calcular_theil_nativo(renta_equiv, f_elev),
        .groups = "drop"
      ) %>%
      mutate(sh_renta = renta_tot / sum(renta_tot))
  
    theil_intra <- sum(tabla_regiones$sh_renta * tabla_regiones$Theil_Interno, na.rm = TRUE)
    theil_inter <- theil_total - theil_intra
    
    return(tibble(
      Año = val_anio,
      Theil_Total = round(theil_total, 4),
      Theil_Intra_CCAA = round(theil_intra, 4),
      Theil_Inter_CCAA = round(theil_inter, 4),
      Pct_Explicado_Brecha_Inter = round((theil_inter / theil_total) * 100, 2)
    ))
  } else {
    return(NULL)
  }
}

tabla_descomposicion_historica <- map_df(anios_estudio, ~ejecutar_descomposicion_anual(.x))

print("--- TABLA EVOLUTIVA DE LA DESCOMPOSICIÓN DE THEIL (CCAA 2016-2021) ---")
print(as.data.frame(tabla_descomposicion_historica))


# ################################################################################
# --- RECONSTRUCCIÓN FINAL HISTÓRICA 
# ################################################################################

# 1. LEER EL ARCHIVO PUENTE 
df_puente <- read_delim(file.path(ruta_base, "Modulo de vidas laborales/MIEF_PERSONA_IDEN.TXT"), 
                           delim=";", col_names=FALSE, show_col_types=FALSE) %>%
  select(ID_SS = X1, ID_HACIENDA_RAW = X8) %>% 
  mutate(
    ID_SS = str_trim(as.character(ID_SS)),
    ID_HACIENDA = as.character(as.integer(str_extract(str_trim(ID_HACIENDA_RAW), "\\d+$")))
  ) %>%
  filter(!is.na(ID_HACIENDA))

# ################################################################################
# --- SELECCIÓN DEL EMPLEO PRINCIPAL POR INDIVIDUO ---
# ################################################################################

ruta_mcvl <- file.path(ruta_base, "Modulo de vidas laborales/MIEF_AFILIAD1.TXT")

data_mcvl <- read_delim(ruta_mcvl, delim=";", col_names=FALSE, show_col_types=FALSE) %>%
  select(ID_SS = X1, GRUPO = X3, CONTRATO = X16) %>%
  mutate(
    ID_SS = str_trim(as.character(ID_SS)),
    contrato_limpio = str_remove(str_trim(CONTRATO), "^0+"),
    primer_digito_c = substr(contrato_limpio, 1, 1),
    
    # 1. Asignamos la estabilidad según metadatos oficiales (1-3 Indefinido | 4-5 Temporal)
    estabilidad_raw = case_when(
      primer_digito_c %in% c("1", "2", "3") ~ "Indefinido",
      primer_digito_c %in% c("4", "5")      ~ "Temporal",
      TRUE                                  ~ "Otros"
    ),
    
    # 2. Asignamos la cualificación según metadatos oficiales (1-3 Alta | 4-11 Media-Baja)
    grupo_limpio = as.numeric(str_extract(GRUPO, "\\d+")),
    cualificacion_raw = case_when(
      grupo_limpio >= 1 & grupo_limpio <= 3  ~ "Alta",
      grupo_limpio >= 4 & grupo_limpio <= 11 ~ "Media-Baja",
      TRUE                                   ~ "Otras"
    )
  ) %>%
  filter(estabilidad_raw != "Otros" & cualificacion_raw != "Otras") %>%
  
group_by(ID_SS) %>%
  summarise(
    estabilidad = if_else(any(estabilidad_raw == "Indefinido"), "Indefinido", "Temporal"),
    cualificacion = if_else(any(cualificacion_raw == "Alta"), "Alta", "Media-Baja"),
    .groups = "drop"
  ) %>%
  mutate(across(c(estabilidad, cualificacion), as.factor))


# EXTRACCIÓN DEL TIPO DE HOGAR
file_21_raw <- list.files(file.path(ruta_base, 2021), pattern = "IDEN", full.names = TRUE)
df_tipohog <- read_lines(file_21_raw) %>%
  tibble(linea = .) %>%
  mutate(
    IDENPER = as.character(as.integer(str_trim(substr(linea, 1, 11)))),
    TIPOHOG = str_extract(linea, "\\d\\.\\d\\.\\d")
  ) %>%
  filter(!is.na(TIPOHOG) & !is.na(IDENPER)) %>%
  select(IDENPER, TIPOHOG)

# ################################################################################
# --- UNIÓN ---
# ################################################################################

base_logit_real_completa <- base_analisis %>%
  filter(anio == 2021) %>%
  mutate(IDENPER = as.character(as.integer(str_trim(IDENPER)))) %>%
  inner_join(df_puente, by = c("IDENPER" = "ID_HACIENDA")) %>%
  inner_join(data_mcvl, by = "ID_SS") %>%
  inner_join(df_tipohog, by = "IDENPER") %>% 
  
  # Ordenamos para priorizar el contrato estable
  arrange(IDENPER, estabilidad) %>%
  distinct(IDENPER, .keep_all = TRUE) %>%
  
  mutate(
    es_pobre = as.numeric(renta_equiv <= quantile(renta_equiv, 0.20, na.rm = TRUE)),
    
    cat_hogar = case_when(
      TIPOHOG == "1.1.1" ~ "1_Unipersonal_Joven",
      TIPOHOG == "1.1.2" ~ "2_Unipersonal_Senior",
      TIPOHOG %in% c("2.1.1", "2.1.2", "2.1.3") ~ "3_Familias_Menores",
      TIPOHOG %in% c("2.2.1", "2.2.2") ~ "4_Adultos_Mayores",
      TRUE ~ "5_Otros_Multicanal"
    ),
    across(c(estabilidad, cualificacion, cat_hogar), as.factor)
  )


# ################################################################################
# --- ANÁLISIS EXPLORATORIO UNIVARIANTE PONDERADO (EJERCICIO 2021) ---
# ################################################################################

# Diseño de encuesta específico para 2021 
diseno_descriptivo_21 <- svydesign(ids = ~1, weights = ~f_elev, data = base_logit_real_completa)

# ==========================================
# A. VARIABLE CONTINUA: RENTA DISPONIBLE EQUIVALENTE
# ==========================================
cat("\n--- ESTADÍSTICOS DESCRIPTIVOS DE LA RENTA EQUIVALENTE PONDERADA ---\n")
media_renta <- svymean(~renta_equiv, diseno_descriptivo_21)
cuartiles_renta <- svyquantile(~renta_equiv, diseno_descriptivo_21, c(0.25, 0.50, 0.75))
print(paste("Media Poblacional:", round(coef(media_renta), 2), "€"))
print("Cuartiles Poblacionales (P25, P50/Mediana, P75):")
print(cuartiles_renta)

# Calculamos la asimetría y curtosis (desviación respecto a la normalidad)
cat("\n--- MOMENTOS DE LA DISTRIBUCIÓN ---\n")
asimetria <- sum(((base_logit_real_completa$renta_equiv - median(base_logit_real_completa$renta_equiv))^3) * base_logit_real_completa$f_elev) / 
  (sum(base_logit_real_completa$f_elev) * sd(base_logit_real_completa$renta_equiv)^3)
cat("Coeficiente de Asimetría de la muestra:", round(asimetria, 3), "\n")

# ==========================================
# B. VARIABLES CATEGÓRICAS: DISTRIBUCIÓN DE LOS PREDICTORES
# ==========================================
cat("\n--- DISTRIBUCIONES DE FRECUENCIAS PONDERADAS (PORCENTAJES POBLACIONALES) ---\n")

# 1. Estabilidad Contractual
cat("\n1. Estabilidad Contractual:\n")
prop_estabilidad <- svymean(~estabilidad, diseno_descriptivo_21)
print(round(prop_estabilidad * 100, 2))

# 2. Cualificación del Puesto
cat("\n2. Cualificación del Puesto:\n")
prop_cualificacion <- svymean(~cualificacion, diseno_descriptivo_21)
print(round(prop_cualificacion * 100, 2))

# 3. Estructura del Hogar
cat("\n3. Estructura del Hogar (5 Categorías):\n")
prop_hogar <- svymean(~cat_hogar, diseno_descriptivo_21)
print(round(prop_hogar * 100, 2))

# ==========================================
# C. VISUALIZACIÓN DESCRIPTIVA
# ==========================================
print(
  ggplot(base_logit_real_completa, aes(x = renta_equiv, weight = f_elev)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_density(color = "darkred", lwd = 1.2) +
    scale_x_continuous(labels = scales::comma, limits = c(0, 75000)) +
    labs(title = "Análisis Univariante: Densidad Poblacional de la Renta Equivalente (2021)",
         x = "Renta Disponible Equivalente (€)", y = "Densidad") +
    theme_minimal()
)

# ################################################################################
# --- ANÁLISIS EXPLORATORIO EVOLUTIVO LONGITUDINAL (2016-2021) ---
# ################################################################################

cat("\n==================================================================\n")
cat("  PARTE A: EVOLUCIÓN DE INDICADORES DE RENTA PONDERADOS (2016-2021)\n")
cat("==================================================================\n")

# Diseñamos el entorno de encuesta para toda la base longitudinal
diseno_longitudinal <- svydesign(ids = ~1, weights = ~f_elev, data = base_analisis)

# ==============================================================================
#   PARTE A: EVOLUCIÓN DE INDICADORES DE RENTA PONDERADOS (2016-2021)
# ==============================================================================
print("--- TABLA EVOLUTIVA DE LA RENTA EQUIVALENTE EN ESPAÑA (CORREGIDA) ---")

# 1. Función para calcular cuantiles ponderados 
fun_cuantil_ponderado <- function(x, w, p = 0.5) {
  ind_validos <- !is.na(x) & !is.na(w)
  x <- x[ind_validos]
  w <- w[ind_validos]
  
  orden <- order(x)
  x_ordenado <- x[orden]
  w_ordenado <- w[orden]
  
  w_acumulado <- cumsum(w_ordenado) / sum(w_ordenado)
  
  # Retorna el valor exacto donde el peso acumulado supera el percentil buscado
  map_dbl(p, function(prob) x_ordenado[which(w_acumulado >= prob)[1]])
}

# 2. Agregación limpia usando la función matemática nativa
resumen_renta_historica <- base_analisis %>%
  group_by(anio) %>%
  summarise(
    Media_Ponderada   = round(weighted.mean(renta_equiv, f_elev, na.rm = TRUE), 2),
    Mediana_Ponderada = round(fun_cuantil_ponderado(renta_equiv, f_elev, 0.50), 2),
    P25_Ponderado     = round(fun_cuantil_ponderado(renta_equiv, f_elev, 0.25), 0),
    P75_Ponderado     = round(fun_cuantil_ponderado(renta_equiv, f_elev, 0.75), 0),
    Muestra_Bruta_N   = n(),
    .groups = "drop"
  )

print(resumen_renta_historica)


cat("\n==================================================================\n")
cat("  PARTE B: EVOLUCIÓN DE PREDICTORES LABORALES (2016-2021)\n")
cat("==================================================================\n")


evidencias_laborales_mcvl <- function(val_anio) {

  file_raw_anio <- list.files(file.path(ruta_base, val_anio), pattern = "IDEN", full.names = TRUE)[1]
  
  df_hogar_anio <- read_lines(file_raw_anio) %>%
    tibble(linea = .) %>%
    mutate(
      IDENPER = as.character(as.integer(str_trim(substr(linea, 1, 11)))),
      TIPOHOG = str_extract(linea, "\\d\\.\\d\\.\\d")
    ) %>%
    filter(!is.na(TIPOHOG) & !is.na(IDENPER))
  
  # Cruzamos la población declarante de ese año con las vidas laborales agregadas por individuo
  base_cruce_anio <- base_analisis %>%
    filter(anio == val_anio) %>%
    mutate(IDENPER = as.character(as.integer(str_trim(IDENPER)))) %>%
    inner_join(df_puente, by = c("IDENPER" = "ID_HACIENDA")) %>%
    inner_join(data_mcvl, by = "ID_SS") %>%
    inner_join(df_hogar_anio, by = "IDENPER")
  
  # Calculamos porcentajes ponderados descriptivos directos
  if(nrow(base_cruce_anio) > 0) {
    diseno_temp <- svydesign(ids = ~1, weights = ~f_elev, data = base_cruce_anio)
    prop_est <- round(svymean(~estabilidad, diseno_temp) * 100, 2)
    prop_cual <- round(svymean(~cualificacion, diseno_temp) * 100, 2)
    
    return(tibble(
      Año = val_anio,
      Muestra_Cruzada = nrow(base_cruce_anio),
      Pct_Indefinidos = prop_est["estabilidadIndefinido"],
      Pct_Temporales = prop_est["estabilidadTemporal"],
      Pct_Cualificacion_Alta = prop_cual["cualificacionAlta"],
      Pct_Cualificacion_Baja = prop_cual["cualificacionMedia-Baja"]
    ))
  } else {
    return(tibble(Año = val_anio, Muestra_Cruzada = 0, Pct_Indefinidos = NA, Pct_Temporales = NA, Pct_Cualificacion_Alta = NA, Pct_Cualificacion_Baja = NA))
  }
}

# Ejecutamos el análisis exploratorio laboral para todos los años de la serie
tabla_exploratoria_longitudinal <- map_df(anios_estudio, ~evidencias_laborales_mcvl(.x))

print("--- TABLA EVOLUTIVA DE LA ESTRUCTURA LABORAL MUESTRAL (2016-2021) ---")
print(tabla_exploratoria_longitudinal)

# ################################################################################
# --- ANÁLISIS EXPLORATORIO UNIVARIANTE AVANZADO
# ################################################################################

# Inicializamos el entorno de encuesta definitivo de 2021 (Muestra limpia individual)
diseno_descriptivo_21 <- svydesign(ids = ~1, weights = ~f_elev, data = base_logit_real_completa)

# ==============================================================================
# TABLA MAESTRA 1: ESTADÍSTICOS DE LA DISTRIBUCIÓN DE LA RENTA EQUIVALENTE
# ==============================================================================
cat("\n==================================================================\n")
cat(" TABLA 1: DESCRIPTIVOS DE LA RENTA DISPONIBLE EQUIVALENTE \n")
cat("==================================================================\n")

# Extraemos los momentos y la dispersión ponderada
media_obj <- svymean(~renta_equiv, diseno_descriptivo_21)
media_val <- as.numeric(coef(media_obj))
media_se  <- as.numeric(SE(media_obj))

# Cuantiles ponderados
q25 <- fun_cuantil_ponderado(base_logit_real_completa$renta_equiv, base_logit_real_completa$f_elev, 0.25)
q50 <- fun_cuantil_ponderado(base_logit_real_completa$renta_equiv, base_logit_real_completa$f_elev, 0.50)
q75 <- fun_cuantil_ponderado(base_logit_real_completa$renta_equiv, base_logit_real_completa$f_elev, 0.75)

# Desviación típica ponderada 
sd_ponderada <- sqrt(sum(base_logit_real_completa$f_elev * (base_logit_real_completa$renta_equiv - media_val)^2) / sum(base_logit_real_completa$f_elev))

tabla_renta_apa <- tibble(
  `Estadístico Descriptivo` = c("Media Aritmética Ponderada", "Error Estándar de la Media (SE)", 
                                "Desviación Típica Ponderada (SD)", "Primer Cuartil (P25)", 
                                "Mediana / Segundo Cuartil (P50)", "Tercer Cuartil (P75)", 
                                "Tamaño de la Muestra Ocupada (N)"),
  `Valor Indicador (€ / Casos)` = c(round(media_val, 2), round(media_se, 2), round(sd_ponderada, 2), 
                                    round(q25, 2), round(q50, 2), round(q75, 2), nrow(base_logit_real_completa))
)

print(as.data.frame(tabla_renta_apa))

# ################################################################################
# --- INFERENCIA DE DISEÑO, CONTRASTES RAO-SCOTT Y DENSIDADES ---
# ################################################################################

# 1. Extracción del puente de sexo para la base descriptiva intermedia
df_sexo_exploratorio <- read_delim(file.path(ruta_base, "Modulo de vidas laborales/MIEF_PERSONA_IDEN.TXT"), 
                                   delim=";", col_names=FALSE, show_col_types=FALSE) %>%
  select(ID_HACIENDA_RAW = X8, COD_SEXO = X3) %>%
  mutate(
    ID_HACIENDA = as.character(as.integer(str_extract(str_trim(ID_HACIENDA_RAW), "\\d+$"))),
    Sexo = case_when(COD_SEXO == "1" ~ "Hombres", COD_SEXO == "2" ~ "Mujeres", TRUE ~ "No_Consta")
  ) %>%
  distinct(ID_HACIENDA, .keep_all = TRUE)

base_descriptiva_completa <- base_logit_real_completa %>%
  inner_join(df_sexo_exploratorio %>% select(ID_HACIENDA, Sexo), by = c("IDENPER" = "ID_HACIENDA"))

diseno_descriptivo_21_ok <- svydesign(ids = ~1, weights = ~f_elev, data = base_descriptiva_completa)


cat("\n==================================================================\n")
cat("  TABLA 2: INFERENCIA POBLACIONAL Y INTERVALOS DE CONFIANZA (95%)  \n")
cat("==================================================================\n")

# Calculamos medias y errores estándar bajo diseño complejo para extraer los IC
obtener_metricas_ic <- function(formula_var, diseno) {
  res_mean <- svymean(formula_var, diseno)
  res_ic   <- confint(res_mean, level = 0.95)
  
  tibble(
    Estimacion = as.numeric(res_mean) * 100,
    Error_Std  = as.numeric(SE(res_mean)) * 100,
    IC_Inf     = as.numeric(res_ic[, 1]) * 100,
    IC_Sup     = as.numeric(res_ic[, 2]) * 100
  )
}

m_est  <- obtener_metricas_ic(~estabilidad, diseno_descriptivo_21_ok)
m_cual <- obtener_metricas_ic(~cualificacion, diseno_descriptivo_21_ok)
m_sex  <- obtener_metricas_ic(~Sexo, diseno_descriptivo_21_ok)
m_hog  <- obtener_metricas_ic(~cat_hogar, diseno_descriptivo_21_ok)

# Unificamos todas las filas métricas 
tabla_inferencial_apa <- tibble(
  `Dimensión Estructural` = c("Estabilidad Contractual", "", "Cualificación del Puesto", "", 
                              "Brecha de Género", "", "Estructura de Convivencia", "", "", "", ""),
  `Categoría Analizada`   = c("Indefinido", "Temporal", "Alta Cualificación", "Media-Baja Cualificación", 
                              "Hombres", "Mujeres", "1_Unipersonal_Joven", "2_Unipersonal_Senior", 
                              "3_Familias_Menores", "4_Adultos_Mayores", "5_Otros_Multicanal"),
  `Porcentaje (%)`        = round(c(m_est$Estimacion, m_cual$Estimacion, m_sex$Estimacion, m_hog$Estimacion), 2),
  `Error Estándar (SE)`   = round(c(m_est$Error_Std, m_cual$Error_Std, m_sex$Error_Std, m_hog$Error_Std), 3),
  `IC Inferior (95%)`     = round(c(m_est$IC_Inf, m_cual$IC_Inf, m_sex$IC_Inf, m_hog$IC_Inf), 2),
  `IC Superior (95%)`     = round(c(m_est$IC_Sup, m_cual$IC_Sup, m_sex$IC_Sup, m_hog$IC_Sup), 2)
)

print(as.data.frame(tabla_inferencial_apa))


cat("\n==================================================================\n")
cat("  TABLA 3: CONTRASTES DE INDEPENDENCIA ASOCIATIVA DE RAO-SCOTT     \n")
cat("==================================================================\n")

# Ejecutamos el test de Rao-Scott 
test_estabilidad   <- svychisq(~estabilidad + es_pobre, diseno_descriptivo_21_ok, statistic = "Chisq")
test_cualificacion <- svychisq(~cualificacion + es_pobre, diseno_descriptivo_21_ok, statistic = "Chisq")
test_hogar         <- svychisq(~cat_hogar + es_pobre, diseno_descriptivo_21_ok, statistic = "Chisq")

tabla_contrastes_apa <- tibble(
  `Asociación del Predictor vs Pobreza` = c("Estabilidad Contractual x Quintil Pobre", 
                                            "Cualificación del Puesto x Quintil Pobre", 
                                            "Estructura del Hogar x Quintil Pobre"),
  `Estadístico F de Rao-Scott`          = c(round(test_estabilidad$statistic, 3), 
                                            round(test_cualificacion$statistic, 3), 
                                            round(test_hogar$statistic, 3)),
  `Valor p (Significación)`             = c(format.pval(test_estabilidad$p.value, eps = 0.001), 
                                            format.pval(test_cualificacion$p.value, eps = 0.001), 
                                            format.pval(test_hogar$p.value, eps = 0.001))
)

print(as.data.frame(tabla_contrastes_apa))


# ==============================================================================
# C. VISUALIZACIÓN BIDIMENSIONAL AVANZADA DE LA RENTA PONDERADA
# ==============================================================================
# Gráfico de densidad acumulada de la renta segmentada por la cualificación ocupacional
print(
  ggplot(base_descriptiva_completa, aes(x = renta_equiv, weight = f_elev, fill = cualificacion)) +
    geom_density(alpha = 0.4, color = "black", lwd = 0.8) +
    scale_x_continuous(labels = scales::comma, limits = c(0, 60000)) +
    scale_fill_manual(values = c("Alta" = "darkgreen", "Media-Baja" = "darkred")) +
    labs(
      title = "Análisis Exploratorio: Densidad de la Renta según Cualificación del Puesto (2021)",
      subtitle = "Población ponderada mediante factores de elevación oficiales del IEF",
      x = "Renta Disponible Equivalente (€)", y = "Densidad Poblacional",
      fill = "Cualificación"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12))
)

# Gráfico de densidad acumulada de la renta segmentada por estructura del hogar
library(ggplot2)
library(ggridges)

# Definimos etiquetas 
etiquetas_limpias <- c(
  "1_Unipersonal_Joven"   = "Unipersonal Joven",
  "2_Unipersonal_Senior"  = "Unipersonal Senior",
  "3_Familias_Menores"    = "Familias con Menores",
  "4_Adultos_Mayores"     = "Adultos Mayores",
  "5_Otros_Multicanal"    = "Otros Hogares"
)

colores_pro <- c(
  "1_Unipersonal_Joven"   = "#2a9d8f", # Turquesa
  "2_Unipersonal_Senior"  = "#264653", # Azul oscuro
  "3_Familias_Menores"    = "#e76f51", # Coral
  "4_Adultos_Mayores"     = "#e9c46a", # Dorado suave
  "5_Otros_Multicanal"    = "#7f8c8d"  # Gris neutro
)


ggplot(base_descriptiva_completa, aes(x = renta_equiv, y = cat_hogar, fill = cat_hogar, weight = f_elev)) +
  geom_density_ridges(alpha = 0.75, color = "white", scale = 1.3, rel_min_height = 0.01, lwd = 0.6) +
  scale_x_continuous(labels = scales::label_comma(big.mark = "."), limits = c(0, 60000)) +
  scale_y_discrete(labels = etiquetas_limpias) +
  scale_fill_manual(values = colores_pro, labels = etiquetas_limpias) +
  
  labs(
    title = "Análisis Exploratorio: Distribución de la Renta según Estructura del Hogar",
    subtitle = "Año 2021",
    x = "Renta Disponible Equivalente (€)",
    y = NULL 
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "gray40", size = 10, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank()
  )

# ################################################################################
# --- MODELO LOGIT 1 (SOCIOLABORAL)  ---
# ################################################################################
if(nrow(base_logit_real_completa) > 0) {
  diseno1 <- svydesign(ids = ~1, weights = ~f_elev, data = base_logit_real_completa)
  
  modelo1 <- svyglm(es_pobre ~ estabilidad + cualificacion + cat_hogar, 
                                 design = diseno1, family = quasibinomial(link = "logit"))
  
  print("\n--- RESULTADOS (ODDS RATIOS) ---")
  print(round(exp(coef(modelo1)), 3))
  
  #  CAPACIDAD PREDICTIVA REAL Y CURVA ROC SUAVIZADA
  base_logit_real_completa$prob <- predict(modelo1, type = "response")
  roc_real <- roc(base_logit_real_completa$es_pobre, base_logit_real_completa$prob)
  plot(roc_real, col = "darkblue", lwd = 4, 
       main = paste("Curva ROC Maestra (5 Categorías de Hogar)\nAUC Final =", round(auc(roc_real), 3)))
  abline(a = 0, b = 1, lty = 2, col = "red")
}



# --- PRUEBA EMPÍRICA DE ROTACIÓN LABORAL ---

mcvl_test <- read_delim(file.path(ruta_base, "Modulo de vidas laborales/MIEF_AFILIAD1.TXT"), 
                        delim=";", col_names=FALSE, show_col_types=FALSE, n_max = 500000) %>%
  select(ID_SS = X1, CONTRATO = X16) %>%
  mutate(
    contrato_limpio = str_remove(str_trim(CONTRATO), "^0+"),
    primer_digito_c = substr(contrato_limpio, 1, 1),
    estabilidad = ifelse(primer_digito_c %in% c("1", "2", "3"), "Indefinido", "Temporal")
  ) %>%
  filter(estabilidad %in% c("Indefinido", "Temporal"))

# Contamos cuántas líneas genera cada trabajador en el archivo según su contrato
tabla_rotacion <- mcvl_test %>%
  group_by(ID_SS, estabilidad) %>%
  summarise(Num_Contratos_Firmados = n(), .groups = "drop") %>%
  group_by(estabilidad) %>%
  summarise(
    Media_Contratos_por_Persona = round(mean(Num_Contratos_Firmados), 2),
    Max_Contratos_un_Solo_Trabajador = max(Num_Contratos_Firmados)
  )

print("--- EVIDENCIA DE ROTACIÓN EN LOS MICRODATOS ---")
print(tabla_rotacion)


length(unique(base_logit_real_completa$IDENPER)) == nrow(base_logit_real_completa)
mean(base_logit_real_completa$es_pobre)


# ################################################################################
# --- AUDITORÍA HISTÓRICA: RECUENTO DE PERSONAS ÚNICAS Y SEXO (2021) ---
# ################################################################################

# 1. Leemos el archivo puente recuperando el Sexo (Columna X3)
df_sexo_puente <- read_delim(file.path(ruta_base, "Modulo de vidas laborales/MIEF_PERSONA_IDEN.TXT"), 
                             delim=";", col_names=FALSE, show_col_types=FALSE) %>%
  select(ID_HACIENDA_RAW = X8, COD_SEXO = X3) %>%
  mutate(
    ID_HACIENDA = as.character(as.integer(str_extract(str_trim(ID_HACIENDA_RAW), "\\d+$"))),
    # Descodificación oficial de la MCVL: 1 = Hombre, 2 = Mujer
    Sexo = case_when(
      COD_SEXO == "1" ~ "Hombres",
      COD_SEXO == "2" ~ "Mujeres",
      TRUE            ~ "No_Consta"
    )
  ) %>%
  distinct(ID_HACIENDA, .keep_all = TRUE)

# 2. Cruzamos el sexo con nuestra base limpia 
base_con_sexo <- base_logit_real_completa %>%
  select(IDENPER, f_elev) %>%
  inner_join(df_sexo_puente, by = c("IDENPER" = "ID_HACIENDA"))

# 3. IMPRIMIMOS LAS CIFRAS
total_unicos <- nrow(base_con_sexo)
tabla_sexo_muestra <- table(base_con_sexo$Sexo)
tabla_sexo_poblacion <- round(prop.table(xtabs(f_elev ~ Sexo, data = base_con_sexo)) * 100, 2)

cat("\n==================================================================\n")
cat("       RECUENTO FINAL DE POBLACIÓN ÚNICA Y COMPOSICIÓN (2021)       \n")
cat("==================================================================\n")
cat("Número total de PERSONAS ÚNICAS en el modelo:", total_unicos, "\n\n")

cat("--- DISTRIBUCIÓN EN LA MUESTRA (Cuentas físicas reales) ---\n")
print(tabla_sexo_muestra)

cat("\n--- DISTRIBUCIÓN POBLACIONAL PONDERADA (Porcentajes oficiales) ---\n")
print(tabla_sexo_poblacion)

# --- AUDITORÍA DE NOMBRES DE VARIABLES REALES ---
colnames(base_analisis)

# ################################################################################
# ---MODELO LOGIT 2 (SOCIOLABORAL AMPLIADO Y TERRITORIAL) Y TABLA COMPARATIVA ---
# ################################################################################

# 1. Recuperación y reconstrucción interna de la variable Sexo
df_sexo_puente_f <- read_delim(file.path(ruta_base, "Modulo de vidas laborales/MIEF_PERSONA_IDEN.TXT"), 
                               delim=";", col_names=FALSE, show_col_types=FALSE) %>%
  select(ID_HACIENDA_RAW = X8, COD_SEXO = X3) %>%
  mutate(
    ID_HACIENDA = as.character(as.integer(str_extract(str_trim(ID_HACIENDA_RAW), "\\d+$"))),
    Sexo = case_when(COD_SEXO == "1" ~ "Hombres", COD_SEXO == "2" ~ "Mujeres", TRUE ~ "No_Consta")
  ) %>%
  distinct(ID_HACIENDA, .keep_all = TRUE)

base_logit_genero <- base_logit_real_completa %>%
  inner_join(df_sexo_puente_f %>% select(ID_HACIENDA, Sexo), by = c("IDENPER" = "ID_HACIENDA")) %>%
  mutate(Sexo = factor(Sexo, levels = c("Hombres", "Mujeres")))

# 2. Extracción Detective de la Región (CCAA) y la Edad (EjnacD) desde el archivo bruto 2021
file_21_raw <- list.files(file.path(ruta_base, 2021), pattern = "IDEN", full.names = TRUE)

df_territorial <- read_lines(file_21_raw) %>%
  tibble(linea = .) %>%
  mutate(
    IDENPER = as.character(as.integer(str_trim(substr(linea, 1, 11)))),
    provincia_cod = substr(linea, 23, 24),
    CCAA = case_when(
      provincia_cod %in% c("04", "11", "14", "18", "21", "23", "29", "41") ~ "Andalucía",
      provincia_cod %in% c("22", "44", "50")                             ~ "Aragón",
      provincia_cod %in% c("33")                                         ~ "Asturias",
      provincia_cod %in% c("07")                                         ~ "Baleares",
      provincia_cod %in% c("35", "38")                                   ~ "Canarias",
      provincia_cod %in% c("39")                                         ~ "Cantabria",
      provincia_cod %in% c("05", "09", "24", "34", "37", "40", "42", "47", "49") ~ "Castilla y León",
      provincia_cod %in% c("02", "13", "16", "19", "45")                 ~ "Castilla-La Mancha",
      provincia_cod %in% c("08", "17", "25", "43")                       ~ "Cataluña",
      provincia_cod %in% c("03", "12", "46")                             ~ "C. Valenciana",
      provincia_cod %in% c("06", "10")                                   ~ "Extremadura",
      provincia_cod %in% c("15", "27", "32", "36")                       ~ "Galicia",
      provincia_cod %in% c("28")                                         ~ "Madrid",
      provincia_cod %in% c("30")                                         ~ "Murcia",
      provincia_cod %in% c("31")                                         ~ "Navarra",
      provincia_cod %in% c("01", "20", "48")                             ~ "País Vasco",
      provincia_cod %in% c("26")                                         ~ "La Rioja",
      TRUE                                                               ~ "Otras/Ceuta/Melilla"
    ),
    anio_nacimiento = as.numeric(str_extract(linea, "\\d{4}")),
    edad_calculada = 2021 - anio_nacimiento
  ) %>%
  filter(!is.na(anio_nacimiento) & !is.na(IDENPER)) %>%
  distinct(IDENPER, .keep_all = TRUE) %>%
  select(IDENPER, CCAA, edad_num = edad_calculada)

# 3. Construcción de la matriz  con todas las variables cruzadas
base_logit_maxima <- base_logit_genero %>%
  inner_join(df_territorial, by = "IDENPER") %>%
  mutate(
    CCAA = factor(CCAA, levels = c("Madrid", "Andalucía", "Aragón", "Asturias", "Baleares", "Canarias", 
                                   "Cantabria", "Castilla y León", "Castilla-La Mancha", "Cataluña", 
                                   "C. Valenciana", "Extremadura", "Galicia", "Murcia", "Navarra", 
                                   "País Vasco", "La Rioja", "Otras/Ceuta/Melilla"))
  )

# 4. Estimación del modelo con todos los controles sociodemográficos
diseno2 <- svydesign(ids = ~1, weights = ~f_elev, data = base_logit_maxima)

modelo2 <- svyglm(es_pobre ~ estabilidad + cualificacion + Sexo + edad_num + CCAA + cat_hogar, 
                           design = diseno2, family = quasibinomial(link = "logit"))

# 5. CAPACIDAD PREDICTIVA 
base_logit_maxima$prob_max <- predict(modelo2, type = "response")
roc_max <- roc(base_logit_maxima$es_pobre, base_logit_maxima$prob_max)

plot(roc_max, col = "darkred", lwd = 4, 
     main = paste("Curva ROC Modelo 2 (Sociodemográfico y Territorial)\nAUC Final =", round(auc(roc_max), 3)))
abline(a = 0, b = 1, lty = 2, col = "grey")

round(exp(coef(modelo1)), 4)
round(exp(confint(modelo1)), 4)
round(coef(modelo1), 4)

cat("\n--- INFORME DE AJUSTE ADICIONAL MAESTRO ---\n")

cat("Criterio de Información de Akaike (AIC):", round(modelo1$aic, 2), "\n")

cat("Deviance Residual del Modelo:", round(modelo1$deviance, 2), "\n")
cat("Grados de Libertad Residuales:", modelo1$df.residual, "\n")

round(AIC(modelo1)[1], 2)

round(exp(coef(modelo2)), 4)
round(exp(confint(modelo2)), 4)
round(coef(modelo2), 4)

cat("\n--- INFORME DE AJUSTE ADICIONAL MAESTRO ---\n")

cat("Criterio de Información de Akaike (AIC):", round(modelo2$aic, 2), "\n")

cat("Deviance Residual del Modelo:", round(modelo2$deviance, 2), "\n")
cat("Grados de Libertad Residuales:", modelo2$df.residual, "\n")

round(AIC(modelo2)[1], 2)

# ==============================================================================
# --- TABLA COMPARATIVA ---
# ==============================================================================
cat("\n==================================================================\n")
cat("      TABLA COMPARATIVA DE DETERMINANTES DE LA VULNERABILIDAD (OE3)    \n")
cat("==================================================================\n")

library(broom)
library(dplyr)

# MODELO 1: Sociolaboral
m1 <- tidy(modelo1) %>%
  mutate(
    OR = round(exp(estimate),4),
    
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    Modelo1 = paste0(OR," ",sig)
  ) %>%
  select(term, Modelo1)


# MODELO 2: Completo
m2 <- tidy(modelo2) %>%
  mutate(
    
    OR = round(exp(estimate),4),
    
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    Modelo2 = paste0(OR," ",sig)
    
  ) %>%
  select(term, Modelo2)


# Unir resultados
tabla_comparativa <- full_join(
  m1,
  m2,
  by="term"
)


# Nombres bonitos para el TFG
tabla_comparativa$term <- recode(
  tabla_comparativa$term,
  
  "(Intercept)"="Intercepto (Riesgo Basal)",
  
  "estabilidadTemporal"=
    "Contrato Temporal (vs Indefinido)",
  
  "cualificacionMedia_Baja"=
    "Cualificación Media-Baja (vs Alta)",
  
  "sexoMujer"=
    "Sexo: Mujer (vs Hombre)",
  
  "edad"=
    "Edad del Trabajador (Años)",
  
  "cat_hogarUnipersonalSenior"=
    "Estructura: Unipersonal Sénior",
  
  "cat_hogarFamiliasConMenores"=
    "Estructura: Familias con Menores",
  
  "cat_hogarAdultosMayores"=
    "Estructura: Adultos Mayores",
  
  "cat_hogarOtrosMulticanal"=
    "Estructura: Otros Multicanal"
)


# Variables que no estén en un modelo
tabla_comparativa[is.na(tabla_comparativa)] <- "No Incluida"


cat("\n--- TABLA COMPARATIVA DE MODELOS ---\n")
print(tabla_comparativa)

cat("\n--- NOTA DE CONTROL: Capacidad Predictiva (AUC) ---\n")
cat("AUC Modelo 1 (Solo Laboral):", round(auc(roc_real), 3), "\n")
cat("AUC Modelo 2 (Sociodemográfico + Regional):", round(auc(roc_max), 3), "\n")
