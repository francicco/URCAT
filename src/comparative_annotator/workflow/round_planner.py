from __future__ import annotations


def choose_next_reference_species(species_order, processed_species, unresolved_by_species):
    candidates = [sp for sp in species_order if sp not in processed_species]
    if not candidates:
        return None

    def keyfunc(sp):
        unresolved = len(unresolved_by_species.get(sp, []))
        seen_bonus = 0 if unresolved > 0 else -1
        return (unresolved, seen_bonus)

    return max(candidates, key=keyfunc)
