import numpy as np

# utility functions for Filters


def get_index_mappings_sample_per_species(managers, sps):
    # get various info from the structures about the center atom species and indexing

    # list of the atom types following the order in managers accross the
    # atomic structures
    types = []
    # list the positions of the begining of each structure in arrays that
    # have one row per atom of species sp
    strides_by_sp = {sp: [0] for sp in sps}
    # count the number of atoms of a particular element sp
    global_counter = {sp: 0 for sp in sps}
    # map position in arrays that have one row per atom of species sp to
    # the position in array that has one row per atom
    indices_by_sp = {sp: [] for sp in sps}
    # map species / structure index / global atom index to the atom index
    # in the structure
    map_by_manager = {sp:[{} for ii in range(len(managers)) ] for sp in sps}
    for i_man in range(len(managers)):
        man = managers[i_man]
        counter = {sp: 0 for sp in sps}
        for i_at, at in enumerate(man):
            types.append(at.atom_type)
            if at.atom_type in sps:
                map_by_manager[sp][i_man][global_counter[at.atom_type]] = i_at
                counter[at.atom_type] += 1
                global_counter[at.atom_type] += 1
            else:
                raise ValueError(
                    'Atom type {} has not been specified in fselect: {}'.format(
                    at.atom_type, sps))
        for sp in sps:
            strides_by_sp[sp].append(counter[sp])

    for sp in sps:
        strides_by_sp[sp] = np.cumsum(strides_by_sp[sp])

    for ii, sp in enumerate(types):
        indices_by_sp[sp].append(ii)

    return strides_by_sp, global_counter, map_by_manager, indices_by_sp


def get_index_mappings_sample(managers):
    # get various info from the structures about the center atom species and indexing

    # list the positions of the begining of each structure in arrays that
    # have one row per atom
    strides = [0]
    # count the number of atoms in managers
    global_counter = 0
    # map the structure index / global atom index to the atom index in the
    # structure
    map_by_manager = [{} for ii in range(len(managers))]
    for i_man in range(len(managers)):
        man = managers[i_man]
        counter = 0
        for i_at, _ in enumerate(man):
            map_by_manager[i_man][global_counter] = i_at
            counter += 1
            global_counter += 1
        strides.append(counter)

    strides = np.cumsum(strides)

    return strides, global_counter, map_by_manager


def convert_selected_global_index2perstructure_index_per_species(managers, selected_ids_by_sp,
                                strides_by_sp, map_by_manager, sps):
    # convert selected center indexing into the rascal format

    # list the atom indices (within their structure) that have been selected
    selected_ids = [[] for ii in range(len(managers))]
    for sp in sps:
        ids = convert_selected_global_index2perstructure_index(
            managers, selected_ids_by_sp[sp], strides_by_sp[sp], map_by_manager[sp])
        for ii, selected_idx in zip(ids, selected_ids):
            selected_idx.extend(ii)
    for ii in range(len(selected_ids)):
        selected_ids[ii] = list(np.sort(selected_ids[ii]))
    return selected_ids


def convert_selected_global_index2perstructure_index(managers, selected_ids_global,
                                                        strides, map_by_manager):
    # convert selected center indexing into the rascal format

    # list the atom indices (within their structure) that have been selected
    selected_ids = [[] for ii in range(len(managers))]
    i_manager = 0
    for idx in selected_ids_global:
        carry_on = True
        while carry_on:
            if idx >= strides[i_manager] and idx < strides[i_manager + 1]:
                selected_ids[i_manager].append(map_by_manager[i_manager][idx])
                carry_on = False
            else:
                i_manager += 1
    for ii in range(len(selected_ids)):
        selected_ids[ii] = list(np.sort(selected_ids[ii]).astype(int))
    return selected_ids
