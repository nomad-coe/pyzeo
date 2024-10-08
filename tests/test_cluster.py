def test_cluster():
    from pyzeo.netstorage import AtomNetwork
    from pyzeo.cluster import pruned_highaccuracy_voronoi_network, \
            get_nearest_largest_diameter_highaccuracy_vornode

    atmnet = AtomNetwork.read_from_CSSR("MgO.cssr", rad_file="MgO.rad")
    vornet, ecs, fcs = atmnet.perform_voronoi_decomposition()
    vornet.analyze_writeto_XYZ('mgo_low', 0.4, atmnet)

    # Pruned high accuracy voronoi network
    ha_vornet = pruned_highaccuracy_voronoi_network(atmnet, 0.7)
    ha_vornet.analyze_writeto_XYZ('mgo_high_prune', 0.4, atmnet)

    # Reduced high accuracy voronoi network
    ha_red_vornet = get_nearest_largest_diameter_highaccuracy_vornode(atmnet) 
    ha_red_vornet.analyze_writeto_XYZ('mgo_high_red', 0.4, atmnet)


