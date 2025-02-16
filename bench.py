class Node:

    def __init__(self, left, right, x, label):
        self.l = left
        self.r = right
        self.x = x
        self.label = label

    def __str__(self):
        return f"Node({self.x}, {self.label})"


def ring_distance(limit, x1, x2):
    return min(
        abs(x2 - x1),
        abs(x2 - x1 + limit),
        abs(x2 - x1 - limit)
    )

def argmedian(x):
    return np.argpartition(x, len(x) // 2)[len(x) // 2]


class CylinderTree:

    def __init__(self, loop_max):
        # All coords are (r, theta).
        # second coord is the loopy one
        self.max = loop_max

    def construct(self, data):
        # The labels are just the full data points.
        idx = np.argsort(data[:, 0])
        data = data[idx]

        # Outer layer is the loopy layer
        self.root = self._construct(data)

    def _construct(self, data):
        # Data is sorted by first coordinate (r), but tree takes second coordinate (theta) first.
        # theta_sort_idx has the indices into data.
        if len(data) == 0:
            return None
        if len(data) == 1:
            return Node(None, None, data[0], data)

        # Split the data array by thetas.
        mid_theta_idx = argmedian(data[:, 1])
        mid_theta = data[mid_theta_idx, 1]

        left_data_sort = []
        right_data_sort = []
        for i, dat in enumerate(data):
            theta = dat[1]
            if i == mid_theta_idx:
                continue
            if theta < mid_theta:
                left_data_sort.append(dat)
            else:
                right_data_sort.append(dat)

        left_data_sort = np.array(left_data_sort)
        right_data_sort = np.array(right_data_sort)
        
        return Node(
            self._construct(left_data_sort),
            self._construct(right_data_sort),
            data[mid_theta_idx], data
        )

    def report_right(self, node, x, y0, y1, ret):
        if node is None:
            return
        if node.x[1] > x:
            self.report_right(node.l, x, y0, y1, ret)
            if node.x[0] > y0 and node.x[0] < y1:
                ret.append(node.x)
            if node.r:
                left = np.searchsorted(node.r.label[:, 0], y0, side='left')
                right = np.searchsorted(node.r.label[:, 0], y1, side='right')
                if left < right:
                    ret.extend(node.r.label[left:right, :])
        else:
            self.report_right(node.r, x, y0, y1, ret)

    def report_left(self, node, x, y0, y1, ret):
        if node is None:
            return
        if node.x[1] < x:
            if node.l:
                left = np.searchsorted(node.l.label[:, 0], y0, side='left')
                right = np.searchsorted(node.l.label[:, 0], y1, side='right')
                if left < right:
                    ret.extend(node.l.label[left:right])
            if node.x[0] > y0 and node.x[0] < y1:
                ret.append(node.x)
            self.report_left(node.r, x, y0, y1, ret)
        else:
            self.report_left(node.l, x, y0, y1, ret)

    def get_range(self, r_range, theta_range):
        r_start, r_end = r_range
        theta_start, theta_end = theta_range
        # Always sweeping from start -> end.
        if theta_start > theta_end:
            # Actually easy case. We just report everything to the right of x_start and left of x_end.
            ret = []
            self.report_right(self.root, theta_start, r_start, r_end, ret)
            self.report_left(self.root, theta_end, r_start, r_end, ret)
            return ret

        node = self.root
        while node:
            if node.x[1] > theta_start and node.x[1] < theta_end:
                ret = []
                self.report_right(node.l, theta_start, r_start, r_end, ret)
                if node.x[0] > r_start and node.x[0] < r_end:
                    ret.append(node.x)
                self.report_left(node.r, theta_end, r_start, r_end, ret)
                return ret
            if node.x[1] < theta_start:
                node = node.r
            elif node.x[1] > theta_end:
                node = node.l
        return []

    def cc_right(self, node, x, y0, y1):
        if node is None:
            return False
        if node.x[1] > x:
            if node.x[0] > y0 and node.x[0] < y1:
                return True
            if node.r:
                left = np.searchsorted(node.r.label[:, 0], y0, side='left')
                right = np.searchsorted(node.r.label[:, 0], y1, side='right')
                if left < right:
                    return True
            return self.cc_right(node.l, x, y0, y1)
        else:
            return self.cc_right(node.r, x, y0, y1)

    def cc_left(self, node, x, y0, y1):
        if node is None:
            return False
        if node.x[1] < x:
            if node.x[0] > y0 and node.x[0] < y1:
                return True
            if node.l:
                left = np.searchsorted(node.l.label[:, 0], y0, side='left')
                right = np.searchsorted(node.l.label[:, 0], y1, side='right')
                if left < right:
                    return True
            return self.cc_left(node.r, x, y0, y1)
        else:
            return self.cc_left(node.l, x, y0, y1)

    def contains_point(self, r_range, theta_range):
        r_start, r_end = r_range
        theta_start, theta_end = theta_range
        # Always sweeping from start -> end.
        if theta_start > theta_end:
            # Actually easy case. We just report everything to the right of x_start and left of x_end.
            return (
                self.cc_right(self.root, theta_start, r_start, r_end)
                or self.cc_left(self.root, theta_end, r_start, r_end)
            )

        node = self.root
        while node:
            if node.x[1] > theta_start and node.x[1] < theta_end:
                if node.x[0] > r_start and node.x[0] < r_end:
                    return True
                return (
                    self.cc_right(node.l, theta_start, r_start, r_end)
                    or self.cc_left(node.r, theta_end, r_start, r_end)
                )
            if node.x[1] < theta_start:
                node = node.r
            elif node.x[1] > theta_end:
                node = node.l
        return False


import time
def test_search_nearest(tree, data, labels, limit, niter=10000):
    numpy_total_time = 0
    rangetree_total_time = 0
    for i in range(niter):
        x_targ = np.random.random()
        t0 = time.time()
        distances = np.min(np.stack([
            np.abs(data - x_targ),
            np.abs(data - x_targ + limit),
            np.abs(data - x_targ - limit)
        ]), axis=0)
        target_label = labels[np.argmin(distances)]
        t1 = time.time()

        binsearch_res = tree.binsearch(x_targ)
        t2 = time.time()
        print(target_label, binsearch_res)
        if target_label != binsearch_res:
            print(x_targ, data[target_label], data[binsearch_res])
            print(ring_distance(limit, x_targ, data[target_label]),
                ring_distance(limit, x_targ, data[binsearch_res]))
        numpy_total_time += t1 - t0
        rangetree_total_time += t2 - t1
    print(numpy_total_time, rangetree_total_time)

def test_range_search(tree, data, labels, limit, niter=1000):
    numpy_total_time = 0
    rangetree_total_time = 0

    idx = np.argsort(data)
    sort_data = data[idx]
    for i in range(niter):
        x_start = np.random.random()
        x_end = np.random.random()
        t0 = time.time()
        if x_start > x_end:
            res = np.concatenate((sort_data[(sort_data > x_start)], sort_data[(sort_data < x_end)]))
        else:
            res = sort_data[(sort_data > x_start) * (sort_data < x_end)]
        t1 = time.time()

        tree_res = np.array([node.x for node in tree.get_range(x_start, x_end)])
        t2 = time.time()
        if len(res) == 0:
            print(len(tree_res))
        else:
            print(np.max(tree_res - res))
        numpy_total_time += t1 - t0
        rangetree_total_time += t2 - t1
    print(numpy_total_time, rangetree_total_time)


def is_in_box(limit, r0, r1, t0, t1, pt):
    if pt[0] < r0 or pt[0] > r1:
        return False
    if t0 < t1:
        return pt[1] > t0 and pt[1] < t1

    return pt[1] > t0 or pt[1] < t1

if __name__ == "__main__":
    import numpy as np
    import ring_rangetree
    np.random.seed(0)
    n_trials = 10
    construct_time = 0
    compute_time = 0
    collide_time = 0
    for k in range(n_trials):
        n = 3000
        data = np.random.random(n)
        labels = np.array(range(n))

        limit = 1

        data2 = np.random.random((n, 2))
        cyl_tree = ring_rangetree.CylinderTree(limit)
        #cyl_tree = CylinderTree(limit)
        t0 = time.time()
        cyl_tree.construct(data2)
        t1 = time.time()
        construct_time += t1 - t0

        for j in range(3000):
            r_start = np.random.random()
            r_end = np.random.random()
            t_start = np.random.random()
            t_end = np.random.random()

            if r_start > r_end:
                tmp = r_end
                r_end = r_start
                r_start = tmp

            t0 = time.time()
            points = cyl_tree.get_range((r_start, r_end), (t_start, t_end))
            t1 = time.time()
            collide = cyl_tree.contains_point((r_start, r_end), (t_start, t_end))
            t2 = time.time()

            compute_time += t1 - t0
            collide_time += t2 - t1

    print("construct:", construct_time)
    print("contain_point:", compute_time)
    print("collide:", collide_time)
