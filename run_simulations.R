rm(list = ls())
pacman::p_load(tidyverse, readxl, lubridate, readr, zoo, cmdstanr, deSolve, 
               wesanderson, foreach, doParallel)

# -----------------------------------------------------------------------------
# DATA DOWNLOAD
# -----------------------------------------------------------------------------

#Tamaño de la población según INEGI
N <- 126014024

if (file.exists("dats_covid.rds") & file.exists("diccionario_covid.rds")){
  message("Reading data from file")
  diccionario.covid <- read_rds("diccionario_covid.rds")
  dats              <- read_rds("dats_covid.rds")
} else {
  message("Downloading data")
  #Descarga de la base principal
  site.covid <- paste0("http://datosabiertos.salud.gob.mx/gobmx/salud",
                       "/datos_abiertos/datos_abiertos_covid19.zip")
  temp <- tempfile()
  download.file(site.covid, temp, method = "curl")
  dats <- read_csv(unz(temp, unzip(temp, list = TRUE)["Name"]))
  unlink(temp)
  
  #Descarga de diccionario de datos para ver el nombre del estado
  temp <- tempfile()
  site.covid.dic <- paste0("http://datosabiertos.salud.gob.mx/gobmx/salud/",
                           "datos_abiertos/diccionario_datos_covid19.zip")
  download.file(site.covid.dic, temp, method = "curl")
  filenames <- unzip(zipfile=temp, list = TRUE)
  fname     <- filenames[which(str_detect(filenames$Name, "Cat.*logo.*")),"Name"]
  unzip(zipfile=temp, files = fname, exdir=".")
  diccionario.covid <- read_excel(fname,
                                  sheet = "Catálogo de ENTIDADES")
  unlink(temp)
  
  write_rds(diccionario.covid, "diccionario_covid.rds")
  write_rds(dats, "dats_covid.rds")
  
  if (dir.exists("diccionario_datos_covid19")){
    unlink("diccionario_datos_covid19", recursive = TRUE)
  }
}

# -----------------------------------------------------------------------------
# DATA CLEANING
# -----------------------------------------------------------------------------

#Pegar para poder juntar
diccionario.covid <- diccionario.covid %>% 
  dplyr::rename(ENTIDAD_RES = CLAVE_ENTIDAD)

#Juntamos nombre del estado
dats <- dats %>% 
  left_join(diccionario.covid, by = "ENTIDAD_RES")

#Cambiamos las cosas a formato fecha
dats <- dats %>% 
  mutate(across(starts_with("FECHA"), ~ ymd(.)))

#Calculamos los muertos
defunciones <- dats %>% 
  filter(!is.na(FECHA_DEF)) %>%
  filter(CLASIFICACION_FINAL %in% c(1,2,3)) %>%
  group_by(FECHA_DEF) %>%
  tally() %>%
  mutate(n = cumsum(n))

#Hospitalizados
hospitalizados <- dats %>% 
  filter(!is.na(FECHA_INGRESO)) %>%
  filter(TIPO_PACIENTE == 2) %>%
  filter(CLASIFICACION_FINAL %in% c(1,2,3)) %>%
  group_by(FECHA_INGRESO) %>%
  tally() %>%
  mutate(n = cumsum(n))

#Casos confirmados
casos <- dats %>%
  filter(!is.na(FECHA_SINTOMAS)) %>%
  filter(CLASIFICACION_FINAL %in% c(1,2,3)) %>%
  group_by(FECHA_SINTOMAS) %>%
  tally() %>%
  mutate(n = cumsum(n))

#UCI
uci <- dats %>% 
  filter(!is.na(FECHA_INGRESO)) %>%
  filter(TIPO_PACIENTE == 2) %>%
  filter(CLASIFICACION_FINAL %in% c(1,2,3)) %>%
  filter(UCI == 1) %>%
  group_by(FECHA_INGRESO) %>%
  tally() %>%
  mutate(n = cumsum(n))

#Juntamos todos y suavizamos
base_julia <- (defunciones %>% 
                 rename(Defunciones = n) %>% 
                 rename(Fecha = FECHA_DEF)) %>%
  left_join( hospitalizados %>% 
               rename(Hospitalizados = n) %>%
               rename(Fecha = FECHA_INGRESO), by = "Fecha") %>%
  left_join( uci %>% 
               rename(UCI = n) %>%
               rename(Fecha = FECHA_INGRESO), by = "Fecha") %>%
  left_join( casos %>% 
               rename(Positivos = n) %>%
               rename(Fecha = FECHA_SINTOMAS), by = "Fecha") %>%
  mutate(across(c(`Defunciones`:`Positivos`), 
                ~rollmean(., 14, align = "left", na.pad = T))) %>%
  mutate(across(c(`Defunciones`:`Positivos`), ~ ./!!N)) %>%
  mutate(t = as.numeric(difftime(Fecha, min(Fecha), units = "days")))

# -----------------------------------------------------------------------------
# MODEL VARIABLES
# -----------------------------------------------------------------------------
base_julia       <- base_julia %>% 
  filter(Fecha < ymd("2021/01/30")) %>% 
  write_excel_csv("base_julia.csv")

base_julia_model <- base_julia %>% 
  filter(t > 0) %>%
  filter(Fecha < ymd("2021/01/30")) 

y0               <- base_julia %>% 
  filter(t == 0) %>% 
  select(Positivos, Hospitalizados, UCI, Defunciones) %>%
  unlist()

exp0 <- base_julia %>% 
  filter(t == 7) %>% 
  select(Positivos) %>%
  unlist()

data_model <- list(
  N            = nrow(base_julia_model),
  t            = base_julia_model$t,
  y0           = c(1.0 - sum(y0), exp0, 0.3*y0["Positivos"], y0["Positivos"], 
                   y0["Hospitalizados"], y0["UCI"],y0["Defunciones"], 
                   y0["Positivos"], y0["Hospitalizados"], y0["UCI"],
                   rep(0, 4)),
  Infected     = base_julia_model$Positivos,
  Hospitalized = base_julia_model$Hospitalizados,
  ICU          = base_julia_model$UCI,
  Death        = base_julia_model$Defunciones,
  Pop          = N,
  hospital_capacity = 7500,
  icu_capacity      = 750,
  q2_low                        = 0.4,
  q2_up                         = 0.6,
  exposed_quarantine_low        = 0.15,
  exposed_quarantine_up         = 0.35,
  t_pre_first                   = 1:100,
  pre_first_quarantine_length   = 100,
  time_first_quarantine_start   = 100,
  time_first_quarantine_end     = 140,
  t_first                       = (100 + 1):140,
  first_quarantine_length       = 40,
  t_pre_second                  = (140 + 1):220,
  pre_second_quarantine_length  = 80,
  time_second_quarantine_start  = 220,
  time_second_quarantine_end    = 260,
  second_quarantine_length      = 40,
  t_second                      = (220 + 1):260,
  t_post_second                 = (260 + 1):nrow(base_julia_model),
  post_second_quarantine_length = length((260 + 1):nrow(base_julia_model)),
  days_on                       = 4,
  days_off                      = 3,
  periodic_quarantine_start     = 100,
  number_of_periods             = floor((nrow(base_julia_model)-100)/(4 + 3)),
  extra_days                    = nrow(base_julia_model) - (100 + 30*7),
  extra_days_on                 = 4,
  extra_days_off                = 1,
  t_pre_periodic                = 1:100,
  t_periodic                    = 1:7,
  t_periodic_on                 = 1:4,
  t_periodic_off                = 1:3,
  t_post_periodic_on            = 1:4,
  t_post_periodic_off            = as.array(1:1)
)

# -----------------------------------------------------------------------------
# MODEL FITTING
# -----------------------------------------------------------------------------
nchains     <- 1
cpath       <-  "/usr/local/opt/llvm/bin/clang++"
cxx_flags   <- "-O3 -march=native"
cpp_options <- list(cxx_flags    = cxx_flags, 
                    cxx          = cpath, 
                    stan_threads = T)
cmd_model <- cmdstan_model("SIR_Model_Chapter_v3.stan", cpp_options = cpp_options)

initial_vals <- function(chain_id = 1) {
 list(
   
  beta_0  = rnorm(1,0.8276086119485562,0.001) %>% abs(), 
  beta_1  = rnorm(1,0.5515727329250524,0.001) %>% abs(), 
  beta_2  = rnorm(1,0.9205242957230675,0.001) %>% abs(), 
  beta_3  = rnorm(1,0.311219365300846,0.001) %>% abs(),
  beta_12 = rnorm(1,0.32561971100893067,0.001) %>% abs(), 
  beta_23 = rnorm(1,0.2015508991016397,0.001) %>% abs(), 
  eta_E   = rnorm(1,0.9960149806661017,0.001) %>% abs(), 
  gamma_0R = rnorm(1,0.7530968545973915,0.001) %>% abs(), 
  gamma_1R = rnorm(1,0.30647950047778616,0.001) %>% abs(), 
  gamma_2R = rnorm(1,0.8292620607852674,0.001) %>% abs(), 
  gamma_3R = rnorm(1,0.16014664285960167,0.001) %>% abs(), 
  theta_3M = rnorm(1,0.46510693286688437,0.001) %>% abs(), 
  coef  = c(rnorm(1,0.0,0.001) %>% abs(),
            rnorm(1,0.009881995293258896,0.001) %>% abs(),
            rnorm(1,0.18586993515558667,0.001) %>% abs(),
            rnorm(1,0.05386425435891228,0.001) %>% abs(),
            rnorm(1,0.01711630633523908,0.001) %>% abs()
            ),
  p_A   = rnorm(1,0.9752388420669302,0.001) %>% abs(), 
  p_I2  = rnorm(1,0.42646170956461865,0.001) %>% abs(), 
  p_I3  = rnorm(1,0.39268337611472104,0.001) %>% abs(), 
  p_M   = rnorm(1,0.8738793655806338,0.001) %>% abs(), 
  #sigma_infected = rnorm(1, 0.1, 0.001) %>% abs(),
  #sigma_hosp = rnorm(1, 0.1, 0.001) %>% abs(),
  #sigma_icu = rnorm(1, 0.1, 0.001) %>% abs(),
  #sigma_dead = rnorm(1, 0.1, 0.001) %>% abs()
  sigma = rnorm(1, 0.1, 0.001) %>% abs()
  )
}

init_ll <- lapply(1:nchains, function(id) initial_vals(chain_id = id))

#GET MLE
m_sample  <- cmd_model$sample(
  data   = data_model,
  seed   = 23789,
  chains = nchains,
  threads_per_chain = 4,
  init = init_ll,
  parallel_chains = nchains,
  iter_warmup = 200,
  iter_sampling = 200,
)

#Get draws and summarized draws
sample_df  <- m_sample$draws(format = "df")
summarized <- m_sample$summary(NULL,~quantile(.x, probs = c(0.025, 0.5, 0.975), na.rm = T))

#Save model 
m_sample$save_output_files(basename = "output_modelo_ajustado_nov_2021_polinomial_chido_vsims")
m_sample$save_object(file = "modelo_ajustado_nov_2021_polinomial_chido_vsims.rds")

# -----------------------------------------------------------------------------
# MODEL CHECK FITTING
# -----------------------------------------------------------------------------

#Plot predictions
infected <- summarized %>%
  filter(str_detect(variable,"Infected_model_1")) %>%
  select(`2.5%`, `50%`, `97.5%`) %>%
  mutate(time = 1:n()) %>%
  mutate(variable = "Infected")

hospitalized <- summarized %>%
  filter(str_detect(variable,"Hospitalized_model_1")) %>%
  select(`2.5%`, `50%`, `97.5%`) %>%
  mutate(time = 1:n()) %>%
  mutate(variable = "Hospitalized")

icu <- summarized %>%
  filter(str_detect(variable,"ICU_model_1")) %>%
  select(`2.5%`, `50%`, `97.5%`) %>%
  mutate(time = 1:n()) %>%
  mutate(variable = "ICU")

death <- summarized %>%
  filter(str_detect(variable,"Death_model_1"))  %>%
  select(`2.5%`, `50%`, `97.5%`) %>%
  mutate(time = 1:n()) %>%
  mutate(variable = "Death")

sims <- infected %>%
  bind_rows(hospitalized) %>%
  bind_rows(icu) %>%
  bind_rows(death) %>%
  mutate(`50%` = if_else(`50%` < 0, 0, `50%`)) %>%
  mutate(Fecha = min(base_julia$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021)

#Reformat cummulative cases
cases <- base_julia %>%
  rename(Death = Defunciones) %>%
  rename(Hospitalized = Hospitalizados) %>%
  rename(ICU = UCI) %>%
  rename(Infected = Positivos) %>%
  rename(time = t) %>%
  pivot_longer(cols = c(`Death`:`Infected`), names_to = "variable") %>%
  filter(year(Fecha) < 2021)

ggplot(sims) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), alpha = 0.1) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable)) +
  geom_point(aes(x = Fecha, y = N*value), size = 0.5, color = "black", data = cases) + 
  geom_point(aes(x = Fecha, y = N*value, color = variable), size = 0.1, data = cases) + 
  facet_wrap(~variable, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  theme_classic() +
  scale_color_manual("Variable", values = wes_palette("Darjeeling1")) +
  scale_fill_manual("Variable", values = wes_palette("Darjeeling1")) +
  labs(
    x = "Date of 2020",
    y = "Number of cases",
    title = "Fitting resulting from the SEIQR model to Mexican population data",
    subtitle = "Probability intervals at 99%"
  ) +
  coord_cartesian(xlim = c(ymd("2020/04/01"),ymd("2020/12/31")))
ggsave("Model_fit_nov_2021.pdf", width = 8, height = 6)

# -----------------------------------------------------------------------------
# MODEL Scenarios
# -----------------------------------------------------------------------------


#MODEL 1
#---------------------------------------------------------------
infected <- summarized %>%
  filter(str_detect(variable,"y_model_1\\[.*,4]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Infected")

hospitalized <- summarized %>%
  filter(str_detect(variable,"y_model_1\\[.*,5]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Hospitalized")

icu <- summarized %>%
  filter(str_detect(variable,"y_model_1\\[.*,6]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "ICU")

death <- summarized %>%
  filter(str_detect(variable,"y_model_1\\[.*,7]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Death")

model1 <- infected %>%
  bind_rows(hospitalized) %>%
  bind_rows(icu) %>%
  bind_rows(death) %>%
  mutate(Measures = "Under current measures")

ggplot(model1) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), alpha = 0.2) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable), size = 1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  scale_color_manual("Variable", values = wes_palette("Darjeeling1")) +
  scale_fill_manual("Variable", values = wes_palette("Darjeeling1")) +
  labs(
    x = "Date of 2020",
    y = "Number of cases",
    title = "Fitting resulting from the SEIQR model to Mexican population data",
    subtitle = "Probability intervals at 99%"
  ) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(ymd("2020/04/01"),ymd("2020/12/31")))
ggsave("Model_as_is_2021.pdf", width = 8, height = 6)

#MODEL 2
#---------------------------------------------------------------
infected <- summarized %>%
  filter(str_detect(variable,"y_model_2\\[.*,4]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Infected")

hospitalized <- summarized %>%
  filter(str_detect(variable,"y_model_2\\[.*,5]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Hospitalized")

icu <- summarized %>%
  filter(str_detect(variable,"y_model_2\\[.*,6]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "ICU")

death <- summarized %>%
  filter(str_detect(variable,"y_model_2\\[.*,7]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Death")

model2 <- infected %>%
  bind_rows(hospitalized) %>%
  bind_rows(icu) %>%
  bind_rows(death) %>%
  mutate(Measures = "Quarantine for 50% symptomatic infections")

ggplot(model2) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), alpha = 0.2) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), size = 1) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), 
              alpha = 0.2, data = model1) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), 
            size = 1, data = model1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  scale_color_manual("Variable", values = wes_palette("Darjeeling1")) +
  scale_fill_manual("Variable", values = wes_palette("Darjeeling1")) +
  labs(
    x = "Date of 2020",
    y = "Number of cases",
    title = "Fitting resulting from the SEIQR model to Mexican population data",
    subtitle = "Probability intervals at 99%"
  ) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(ymd("2020/04/01"),ymd("2020/12/31")))
ggsave("Model_q2_2021.pdf", width = 8, height = 6)

#MODEL 3
#---------------------------------------------------------------
infected <- summarized %>%
  filter(str_detect(variable,"y_model_3\\[.*,4]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Infected")

hospitalized <- summarized %>%
  filter(str_detect(variable,"y_model_3\\[.*,5]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Hospitalized")

icu <- summarized %>%
  filter(str_detect(variable,"y_model_3\\[.*,6]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "ICU")

death <- summarized %>%
  filter(str_detect(variable,"y_model_3\\[.*,7]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Death")

model3 <- infected %>%
  bind_rows(hospitalized) %>%
  bind_rows(icu) %>%
  bind_rows(death) %>%
  mutate(Measures = "Quarantine for susceptibles")

ggplot(model3) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), alpha = 0.2) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), size = 1) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), 
              alpha = 0.2, data = model1) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), 
            size = 1, data = model1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  scale_color_manual("Variable", values = wes_palette("Darjeeling1")) +
  scale_fill_manual("Variable", values = wes_palette("Darjeeling1")) +
  labs(
    x = "Date of 2020",
    y = "Number of cases",
    title = "Fitting resulting from the SEIQR model to Mexican population data",
    subtitle = "Probability intervals at 99%"
  ) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(ymd("2020/04/01"),ymd("2020/12/31")))
ggsave("Model_q2_2021.pdf", width = 8, height = 6)


#MODEL 5
#---------------------------------------------------------------
infected <- summarized %>%
  filter(str_detect(variable,"y_model_5\\[.*,4]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Infected")

hospitalized <- summarized %>%
  filter(str_detect(variable,"y_model_5\\[.*,5]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Hospitalized")

icu <- summarized %>%
  filter(str_detect(variable,"y_model_5\\[.*,6]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "ICU")

death <- summarized %>%
  filter(str_detect(variable,"y_model_5\\[.*,7]")) %>%
  mutate(time = 1:n()) %>%
  mutate(Fecha = min(cases$Fecha) + days(time)) %>%
  filter(year(Fecha) < 2021) %>%
  mutate(variable = "Death")

model5 <- infected %>%
  bind_rows(hospitalized) %>%
  bind_rows(icu) %>%
  bind_rows(death) %>%
  mutate(Measures = "Model with saturation")

hlinepos <- data.frame(
  yline = c(0,7500,750,0),
  variable = c("Death","Hospitalized","ICU","Infected")
)

ggplot(model5) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), alpha = 0.2) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), size = 1) +
  geom_ribbon(aes(x = Fecha, ymin = N*`2.5%`, ymax = N*`97.5%`, fill = variable), 
              alpha = 0.2, data = model1) +
  geom_line(aes(x = Fecha, y = N*`50%`, color = variable, linetype = Measures), 
            size = 1, data = model1) +
  theme_bw() +
  geom_hline(aes(yintercept = yline), data = hlinepos, linetype = "dotdash") +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  scale_color_manual("Variable", values = wes_palette("Darjeeling1")) +
  scale_fill_manual("Variable", values = wes_palette("Darjeeling1")) +
  labs(
    x = "Date of 2020",
    y = "Number of cases",
    title = "Fitting resulting from the SEIQR model to Mexican population data",
    subtitle = "Probability intervals at 99%"
  ) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(ymd("2020/04/01"),ymd("2020/12/31"))) 
ggsave("Model_sat_2021.pdf", width = 8, height = 6)



