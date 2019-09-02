import json
import numpy as np


def adapt_structure(cell, positions, numbers, pbc):
    cell = np.array(cell.T, order='F')
    positions = np.array(positions.T, order='F')
    numbers = numbers.reshape(-1, 1)
    pbc = pbc.reshape(3, 1)
    return dict(cell=cell, positions=positions, atom_types=numbers, pbc=pbc)


def dump_json_frame(fn, frames):
    from collections import Iterable
    if not isinstance(frames, Iterable):
        frames = [frames]
    data = dict()
    for ii, frame in enumerate(frames):
        data[ii] = dict(positions=frame.get_positions().tolist(),
                        cell=frame.get_cell().tolist(),
                        numbers=frame.get_atomic_numbers().tolist(),
                        pbc=frame.get_pbc().tolist())

    data['ids'] = np.arange(len(frames)).tolist()
    data['nextid'] = 2
    with open(fn, 'w') as f:
        json.dump(data, f, indent=2, separators=(',', ': '))


def load_json_frame(fn):
    with open(fn, 'r') as f:
        data = json.load(f)
    ids = data['ids']
    structure = {key: np.array(val)
                 for idx in ids for key, val in data[str(idx)].items()}

    return adapt_structure(**structure)


class BoxList(object):
    def __init__(self, max_cutoff, cell, pbc, centers):
        # Compute reciprocal lattice vectors.
        b1_c, b2_c, b3_c = np.linalg.pinv(cell).T

        # Compute distances of cell faces (height between 2
        # consecutive faces [010]
        l1 = np.linalg.norm(b1_c)
        l2 = np.linalg.norm(b2_c)
        l3 = np.linalg.norm(b3_c)
        face_dist_c = np.array([1 / l1 if l1 > 0 else 1,
                                1 / l2 if l2 > 0 else 1,
                                1 / l3 if l3 > 0 else 1])

        # We use a minimum bin size of 3 A
        self.bin_size = max_cutoff
        # Compute number of bins such that a sphere of radius
        # cutoff fit into eight neighboring bins.
        self.nbins_c = np.maximum(
            (face_dist_c / self.bin_size).astype(int), [1, 1, 1])
        self.nbins = np.prod(self.nbins_c)
        # Compute over how many bins we need to loop in the
        # neighbor list search.
        self.neigh_search = np.ceil(
            self.bin_size * self.nbins_c / face_dist_c).astype(int)
        self.bin2icenters = [[] for bin_idx in range(self.nbins)]
        scaled_positions_ic = np.linalg.solve(cell.T, centers.T).T
        self.h_sizes = np.linalg.norm(cell, axis=1)
        self.part2bin = {}
        for icenter in range(len(centers)):
            bin_index_ic = np.floor(
                scaled_positions_ic[icenter]*self.nbins_c).astype(int)
            bin_id = self.cell2lin(bin_index_ic)
            self.bin2icenters[bin_id].append(icenter)
            self.part2bin[icenter] = bin_id
            self.list = []
        # print(self.nbins)
        for bin_id in range(self.nbins):
            self.list.append(Box(
                bin_id, self.nbins_c, self.neigh_search,
                self.bin2icenters[bin_id], pbc, self))

    def cell2lin(self, ids):
        return int(ids[0] + self.nbins_c[0] * (ids[1] + self.nbins_c[1]
                    * ids[2]))

    def iter_box(self):
        for bin_id in range(self.nbins):
            yield self.list[bin_id]

    def __getitem__(self, bin_id):
        return self.list[bin_id]


class Box(object):
    def __init__(self, lin_pos, nbins_c, neigh_search, icenters, pbc, boxlist):
        self.nbins_c = nbins_c
        self.neigh_search = neigh_search
        self.icenters = icenters
        self.pbc = pbc
        self.lin_pos = lin_pos
        self.mult_pos = self.lin2cell(lin_pos)
        self.boxlist = boxlist
        self.search_idx = []

        for ii, p in enumerate(self.pbc):

            if 0 == self.mult_pos[ii] and p is False:
                self.search_idx.append(
                    [self.mult_pos[ii]+jj for jj in range(
                                                    self.neigh_search[ii]+1)])
            elif self.nbins_c[ii]-1 == self.mult_pos[ii] and p is False:
                self.search_idx.append(
                    [self.mult_pos[ii]+jj for jj in range(
                                                -self.neigh_search[ii], 0+1)])
            else:
                self.search_idx.append([self.mult_pos[ii]+jj for jj in
                                        range(-self.neigh_search[ii],
                                        self.neigh_search[ii]+1)])

        self.neighbour_bin_index, self.neighbour_bin_shift = [], []
        for ii in self.search_idx[0]:
            for jj in self.search_idx[1]:
                for kk in self.search_idx[2]:
                    box_shift, box_pos = np.divmod([ii, jj, kk], self.nbins_c)
                    neigh_box_idx = self.cell2lin(box_pos)
                    self.neighbour_bin_index.append(neigh_box_idx)
                    self.neighbour_bin_shift.append(box_shift)

    def cell2lin(self, ids):
        return int(ids[0] + self.nbins_c[0] * (ids[1] + self.nbins_c[1]
                    * ids[2]))

    def lin2cell(self, lin_ids):
        fac = 1
        cell_pos = np.array([0, 0, 0])
        for ii in range(3):
            cell_pos[ii] = lin_ids/fac % self.nbins_c[ii]
            fac *= self.nbins_c[ii]
        return cell_pos

    def iter_neigh_box(self):
        from copy import deepcopy
        for ii in self.search_idx[0]:
            for jj in self.search_idx[1]:
                for kk in self.search_idx[2]:
                    box_shift, box_pos = np.divmod([ii, jj, kk], self.nbins_c)
                    neigh_box_idx = self.cell2lin(box_pos)
                    jcenters = deepcopy(self.boxlist[neigh_box_idx].icenters)
                    for jneigh in jcenters:
                        yield jneigh, deepcopy(box_shift)
