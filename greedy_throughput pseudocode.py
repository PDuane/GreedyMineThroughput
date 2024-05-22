dropped_nodes = [starting_node]
for N nodes to drop
    for p1 = each non-obstructed pixel
        best_node = previously placed node with best data rate to this pixel
        if throughput between best_node and p1 > throughput 
            for p2 = each non-obstructed pixel
                coverage_map(p2) = max(coverage_map(p1), throughput between p2 and p1)
        
            if sum(coverage_map > best_coverage)
                best_place = p1
    add best_place to dropped_nodes