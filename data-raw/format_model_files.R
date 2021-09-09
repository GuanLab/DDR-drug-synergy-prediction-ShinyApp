
all_model_path = list.files(path = here::here("inst/models/"),
                            full.names=T,
                            pattern = "aoc")

purrr::walk(all_model_path, function(x) {
  
  model = lightgbm::lgb.load(x)
  
  rds_file = gsub("txt", "RDS", x)
  lightgbm::saveRDS.lgb.Booster(model, rds_file)
})


all_model_path = list.files(path = here::here("inst/models/"),
                            full.names=T,
                            pattern = "bliss")

purrr::walk(all_model_path, function(x) {
  
  model = lightgbm::lgb.load(x)
  
  rds_file = gsub("txt", "RDS", x)
  lightgbm::saveRDS.lgb.Booster(model, rds_file)
})
