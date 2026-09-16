# =============================================================================
# workflow_constants.R
# Constant lookup values used by cleaning functions.
# =============================================================================

trusted_resources <- c(
  "BirdLife Australia, Birdata",
  "NSW BioNet Atlas",
  "Victorian Biodiversity Atlas",
  "SA Fauna",
  "Tasmanian Natural Values Atlas",
  "eBird Australia",
  "iNaturalist Australia"
)

universal_protocol_exclusions <- c(
  "Traveling", "Historical", "Roadkill", "Dead", "Sign", "Scat",
  "Burrow", "Tracks", "Nest", "Eggs", "unknown", "pelagic",
  "Wetland", "Waterbird", "other", "Wetlands", "Translocation"
)

tas_endemics <- c(
  "Acanthiza ewingii", "Acanthornis magna", "Anthochaera paradoxa",
  "Melanodryas vittata", "Melithreptus affinis",
  "Melithreptus validirostris", "Nesoptilotis flavicollis",
  "Pardalotus quadragintus", "Platycercus caledonicus",
  "Sericornis humilis", "Strepera fuliginosa"
)

mainland_obligates <- c(
  "Acanthagenys rufogularis", "Acanthiza lineata", "Acanthiza nana",
  "Acanthiza uropygialis", "Alectura lathami", "Aphelocephala leucopsis",
  "Aprosmictus erythropterus", "Artamus cinereus", "Artamus leucorynchus",
  "Barnardius zonarius", "Caligavis chrysops", "Callocephalon fimbriatum",
  "Calyptorhynchus banksii", "Calyptorhynchus lathami",
  "Centropus phasianinus", "Chlamydera maculata", "Cincloramphus mathewsi",
  "Cinclosoma castanotum", "Climacteris erythrops", "Climacteris picumnus",
  "Coracina maxima", "Coracina papuensis", "Corcorax melanorhamphos",
  "Cormobates leucophaea", "Cracticus nigrogularis",
  "Daphoenositta chrysoptera", "Dicrurus bracteatus",
  "Drymodes brunneopygia", "Edolisoma tenuirostre", "Entomyzon cyanotis",
  "Eopsaltria australis", "Gerygone fusca", "Gerygone levigaster",
  "Gerygone mouki", "Gerygone olivacea", "Grantiella picta",
  "Hylacola cauta", "Hylacola pyrrhopygia", "Lichenostomus cratitius",
  "Lichenostomus melanops", "Malurus lamberti", "Malurus melanocephalus",
  "Malurus splendens", "Manorina flavigula", "Melanodryas cucullata",
  "Melithreptus albogularis", "Melithreptus gularis", "Merops ornatus",
  "Microeca fascinans", "Myiagra inquieta", "Myiagra rubecula",
  "Myzomela sanguinolenta", "Northiella haematogaster", "Oreoica gutturalis",
  "Origma solitaria", "Oriolus sagittatus", "Pachycephala inornata",
  "Pachycephala rufogularis", "Petroica goodenovii",
  "Philemon citreogularis", "Philemon corniculatus",
  "Plectorhyncha lanceolata", "Polytelis anthopeplus",
  "Polytelis swainsonii", "Pomatostomus superciliosus",
  "Pomatostomus temporalis", "Psephotus haematonotus", "Ptilotula fusca",
  "Ptilotula ornata", "Ptilotula penicillata", "Purnella albifrons",
  "Pycnonotus jocosus", "Pyrrholaemus sagittatus",
  "Scythrops novaehollandiae", "Smicrornis brevirostris",
  "Stagonopleura guttata", "Stizoptera bichenovii", "Struthidea cinerea",
  "Sugomel niger", "Taeniopygia guttata", "Todiramphus pyrrhopygius"
)

# Pass the complete rule set explicitly to the cleaning function. Keeping this
# as one object makes the scientific lookup values visible without making the
# function depend on loose objects in the global workspace.
cleaning_rules <- list(
  trusted_resources = trusted_resources,
  protocol_exclusions = universal_protocol_exclusions,
  tasmanian_endemics = tas_endemics,
  mainland_obligates = mainland_obligates
)
