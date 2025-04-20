library(rgbif)

# Download a dataset
options(gbif_user = "", 
        gbif_pwd = "", 
        gbif_email = "")

#########################################################################
# Reference
# https://pacificwild.org/pacific-salmon-species-spotlight/
# https://www.pac.dfo-mpo.gc.ca/fm-gp/salmon-saumon/facts-infos-eng.html
#########################################################################

# Chinook Salmon
name_backbone("Oncorhynchus tshawytscha")
occ_search(scientificName = "Oncorhynchus tshawytscha",country = "CA")
# occ_download(pred("taxonKey", 5204024), format = "SIMPLE_CSV")
gbif_download <- occ_download(pred("taxonKey", 5204024),format = "SIMPLE_CSV")
occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Chum Salmon
name_backbone("Oncorhynchus keta")
gbif_download <- occ_download(pred("taxonKey", 5204014),
                              pred_gt("year", 2009),
                              pred("country","CA"),
                              format = "SIMPLE_CSV")
occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Pink Salmon
name_backbone("Oncorhynchus gorbuscha")
gbif_download <- occ_download(pred("taxonKey", 5204037),
                              pred_gt("year", 2009),
                              pred("country","CA"),
                              format = "SIMPLE_CSV")
occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Coho Salmon
name_backbone("Oncorhynchus kisutch")
gbif_download <- occ_download(pred("taxonKey", 5204034),
                              pred_gt("year", 2009),
                              pred("country","CA"),
                              format = "SIMPLE_CSV")
occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Sockeye Salmon
name_backbone("Oncorhynchus nerka")
gbif_download <- occ_download(pred("taxonKey", 5204039),
                              pred_gt("year", 2009),
                              pred("country","CA"),
                              format = "SIMPLE_CSV")
occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Canada Goose(Branta canadensis)
name_backbone("Branta canadensis")
gbif_download <- occ_download(pred("taxonKey", 5232437),
                              pred_or(pred("year", 2010), pred("year", 2024)),
                              pred("country","CA"),
                              format = "SIMPLE_CSV")

occ_download_wait(gbif_download)
d <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Load datasets
chinook_salmon <- read.csv("data/chinook_salmon.csv", sep="\t")
pink_salmon <- read.csv("data/pink_salmon.csv", sep="\t")
coho_salmon <- read.csv("data/coho_salmon.csv", sep="\t")
chum_salmon <- read.csv("data/chum_salmon.csv", sep="\t")
sockeye_salmon <- read.csv("data/sockeye_salmon.csv", sep="\t")
salmon <- rbind(
  chinook_salmon,
  pink_salmon,
  coho_salmon,
  chum_salmon,
  sockeye_salmon
)
salmon <- salmon[!is.na(salmon$decimalLongitude) & !is.na(salmon$decimalLatitude), ]
salmon <- salmon[salmon$year >=2010,]

write.csv(salmon,
          file   = "data/salmon.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
