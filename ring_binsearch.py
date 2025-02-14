
class Node:

    def __init__(self, left, right, x, label):
        self.l = left
        self.r = right
        self.x = x
        self.label = label

    def __str__(self):
        return f"Node({self.x}, {self.label})"

class RingTree:

    def __init__(self, max_x):
        self.max = max_x
        self.root = None

    def construct(self, data, labels):
        idx = np.argsort(data)
        data = data[idx]
        labels = labels[idx]
        self.root = self._construct(data, labels)

    def _construct(self, data, labels):
        """Returns a root and does so recursively."""
        if len(data) == 0:
            return None
        if len(data) == 1:
            return Node(None, None, data[0], labels[0])
        
        mid_idx = (len(data) - 1) // 2
        return Node(
            self._construct(data[:mid_idx], labels[:mid_idx]),
            self._construct(data[mid_idx+1:], labels[mid_idx+1:]),
            data[mid_idx], labels[mid_idx]
        )

    def binsearch(self, x):
        cand1 = self.find_close_right(self.root, x)
        cand2 = self.find_close_left(self.root, x)

        if cand1 is None:
            cand1 = self.root
            while cand1.l:
                cand1 = cand1.l
            dist1 = cand1.x - x + self.max
        else:
            dist1 = cand1.x - x

        if cand2 is None:
            cand2 = self.root
            while cand2.r:
                cand2 = cand2.r
            dist2 = x - cand2.x + self.max
        else:
            dist2 = x - cand2.x

        if dist1 < dist2:
            return cand1.label
        return cand2.label

    def find_close_right(self, node, x):
        # Find the closest node that is greater or equal to x.
        if node.x == x:
            return node
        if node.x > x:
            if node.l is None:
                return node
            cand = self.find_close_right(node.l, x)
            if cand is None:
                return node
            if cand.x < node.x:
                return cand
            return node

        if node.r is None:
            return None     # Nothing bigger than x.
        return self.find_close_right(node.r, x)

    def find_close_left(self, node, x):
        # Find the closest node that is less or equal to x.
        if node.x == x:
            return node
        if node.x < x:
            if node.r is None:
                return node
            cand = self.find_close_left(node.r, x)
            if cand is None:
                return node
            if cand.x > node.x:
                return cand
            return node

        if node.l is None:
            return None     # Nothing less than x.
        return self.find_close_left(node.l, x)
    

    def report_right(self, node, x, ret):
        if node is None:
            return
        if node.x > x:
            self.report_right(node.l, x, ret)
            ret.append(node)
        self.report_right(node.r, x, ret)
    
    def report_left(self, node, x, ret):
        if node is None:
            return
        self.report_left(node.l, x, ret)
        if node.x < x:
            ret.append(node)
            self.report_left(node.r, x, ret)

    def get_range(self, x_start, x_end):
        # Always sweeping from x_start -> x_end.
        if x_start > x_end:
            # Actually easy case. We just report everything to the right of x_start and left of x_end.
            ret = []
            self.report_right(self.root, x_start, ret)
            self.report_left(self.root, x_end, ret)
            return ret

        node = self.root
        while node:
            if node.x > x_start and node.x < x_end:
                ret = []
                self.report_right(node.l, x_start, ret)
                ret.append(node)
                self.report_left(node.r, x_end, ret)
                return ret
            if node.x < x_start:
                node = node.r
            elif node.x > x_end:
                node = node.l
        return []

    def traverse(self):
        ret = []
        self._traverse(self.root, ret)
        return ret

    def _traverse(self, node, ret):
        if node is None:
            return ret
        self._traverse(node.l, ret)
        ret.append(node.x)
        self._traverse(node.r, ret)

def ring_distance(limit, x1, x2):
    return min(
        abs(x2 - x1),
        abs(x2 - x1 + limit),
        abs(x2 - x1 - limit)
    )


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

if __name__ == "__main__":
    import numpy as np
    np.random.seed(0)
    n = 100
    data = np.random.random(n)
    labels = np.array(range(n))

    limit = 1

    tree = RingTree(limit)
    tree.construct(data, labels)
    res = tree.traverse()

    #test_search_nearest(tree, data, labels, limit, niter=10000)
    test_range_search(tree, data, labels, limit, niter=1000)

