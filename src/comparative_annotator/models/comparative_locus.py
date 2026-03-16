class ComparativeLocus:

    def __init__(self, locus_id, seed_species, seed_transcript):
        self.locus_id = locus_id
        self.seed_species = seed_species
        self.seed_transcript = seed_transcript

        self.primary = {}       # species -> best locus
        self.alternatives = {}  # species -> other candidates

    def set_primary(self, species, locus_id):
        self.primary[species] = locus_id

    def set_alternatives(self, species, locus_ids):
        self.alternatives[species] = locus_ids
