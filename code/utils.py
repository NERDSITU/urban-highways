import warnings
warnings.filterwarnings("ignore")

import sys
import argparse
import random
from time import time
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import LineString, Point, MultiPoint, MultiLineString, Polygon
from tqdm import tqdm
from geopy.distance import distance, geodesic, great_circle
import osmnx as ox
import networkx as nx
import scipy.stats
from scipy.stats import ks_2samp
import sklearn
import igraph
from igraph import Graph

from random import choice
from bisect import bisect_left
import copy

from functools import partial
import pyproj
from pyproj import Geod
from pyproj.crs import ProjectedCRS
from pyproj.crs.coordinate_operation import AzimuthalEquidistantConversion
from shapely.ops import transform

######### FUNCTIONS ######### 

def check_configs(config):
    '''helper function to check that config file is correctly filled out'''
    beeline = config["dist"] 
    density = config["density"]
    rand_layout = config["rand_layout"]
    seed = config["seed"]
    cbsacode = config["cbsa"]
    
    assert beeline in ["walking", "beeline"], "Error: DIST must be one of: 'beeline', 'walking'! Please correct in ./code/config.yml"
    assert density in [0, 1], "Error: DENSITY must be one of: 0, 1! Please correct in ./code/config.yml"
    assert rand_layout in [0, 1], "Error: RAND_LAYOUT must be one of: 'beeline', 'walking'! Please correct in ./code/config.yml"
    assert type(seed)==int, "Error: SEED must be an integer! Please correct in ./code/config.yml"
    assert cbsacode in [12420, 16980, 33100, 38300, 38900], "Error: CBSA must be one of: 12420, 16980, 33100, 38300, 38900! Please correct in ./code/config.yml"
    
    print("Config file filled out correctly, proceeding")
    
    return None


# line for edges
def create_polyline(r):
    return LineString([[r["lon_home1"],r["lat_home1"]],[r["lon_home2"],r["lat_home2"]]])

# point for nodes
def create_point(r):
    return Point([r["lon_home"],r["lat_home"]])

def geodesic_point_buffer(lon, lat, meters):
    proj_crs = ProjectedCRS(conversion = AzimuthalEquidistantConversion(lat, lon))
    proj_wgs84 = pyproj.Proj('EPSG:4326')
    Trans = pyproj.Transformer.from_proj(proj_crs, proj_wgs84, always_xy=True).transform
    return transform(Trans, Point(0, 0).buffer(meters))

def plot_city(nodes_gdf, edges_gdf):
    """
    plot individual level mutual followership network inside metropolitan areas
    """
    # map of social connections -- geolocated nodes and edges
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    nodes_gdf.plot(markersize=3, color='k', ax=ax)
    edges_gdf.plot(lw=0.1, color='grey', alpha=0.2, ax=ax)
    ax.axis('off')
    plt.tight_layout()
    plt.show()
    return (nodes_gdf, edges_gdf)

def collate_street_types(highway, motorway_osm_labels):
    """
    Aggregate similar OSM street types
    """
    if highway in motorway_osm_labels:
        return 'motorway'
    elif highway == 'tertiary':
        return 'secondary'
    else:
        return highway
    
def ci(a):
    if len(a) == 1:
        return 0
    sd = scipy.stats.sem(a)
    avg = np.mean(a)
    if sd == 0:
        return 0
    return scipy.stats.t.interval(0.95, len(a)-1, loc=avg, scale=sd)

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def convert_nxgraph_to_igraph(G):
    """
    Converts networkx to igraph, setting the networkx ids as the name in igraph node attributes
    @param G: a networkx graph
    @return G_igraph: an igraph graph
    """
    G_igraph = Graph.from_networkx(G)
    for v in G_igraph.vs:
        v['name'] = v.attributes()['_nx_name']
    return G_igraph

def get_street_network(
        polygon, 
        customfilter,
        cbsacode,
        dataset,
        savetofolder=None,
    ):

    # get OSM street network
    print(f"Applying custom filter: {customfilter}")
    G = ox.graph.graph_from_polygon(
        polygon=polygon,
        simplify=False, 
        retain_all=True, 
        truncate_by_edge=True, 
        clean_periphery=False, 
        custom_filter=customfilter
        )
    print(f"OSM data downloaded for {cbsacode}")

    G_simp = ox.simplify_graph(
            G = G, # the directed, not yet projected graph 
            edge_attrs_differ=["highway"], 
            strict=True, 
            remove_rings=False
            )
    print("OSM data simplified")

    if savetofolder:
        print("Saving G_simplified...")
        filename = f'{cbsacode}_G_simplified_{dataset}.graphml'
        ox.save_graphml(
            G = G_simp, 
            filepath = savetofolder + filename
        )

    G_simp = ox.project_graph(G_simp, to_crs='epsg:4326')
    G_simp = G_simp.to_undirected()

    # turn the network into a geodataframe
    gdf_nodes_city, gdf_edges_city = ox.graph_to_gdfs(G_simp)
    # gdf_edges_city['saved_geom'] = gdf_edges_city.geometry # ANVY: drop this?
    
    # saving the street network and its lists of nodes and edges
    if savetofolder:

        print("Saving gdfs...")
        
        gdf_nodes_city.to_file(
            savetofolder + f'{cbsacode}_streetnodes_{dataset}.gpkg', index = False
            )
        gdf_edges_city[["geometry"]].to_file(
            savetofolder + f'{cbsacode}_streetedges_{dataset}.gpkg', index = False
            )
        
        # ANVY: drop all pickles below?
        with open(savetofolder + f'{cbsacode}_streetedges_{dataset}.pickle', 'wb') as handle:
            pickle.dump(gdf_edges_city, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(savetofolder + f'{cbsacode}_streetnodes_{dataset}.pickle', 'wb') as handle:
            pickle.dump(gdf_nodes_city, handle, protocol=pickle.HIGHEST_PROTOCOL)    
        with open(savetofolder + f'{cbsacode}_streetnetwork_{dataset}.pickle', 'wb') as handle:
            pickle.dump(G, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return G, gdf_nodes_city, gdf_edges_city

def get_city_geodataframes(nodes_gpkg, edges_gpkg, cbsacode, city_boundaries=None):
    """
    @param nodes_gpkg: gpkg with all social network nodes
    @param edges_gpkg: gpkg with all social network edges
    @param cbsacode: the cbsacode of the area of interest
    @param city_boundaries: polygon determining city boundaries
    @return gdf_nodes, gdf_edges: the dataframes with the info about nodes and edges
    """

    gdf_nodes = gpd.read_file(nodes_gpkg)
    gdf_edges = gpd.read_file(edges_gpkg)

    gdf_nodes = gdf_nodes[gdf_nodes["cbsacode"]==cbsacode].copy().reset_index(drop=True)
    gdf_edges = gdf_edges[gdf_edges["cbsacode"]==cbsacode].copy().reset_index(drop=True)

    # limiting nodes and edges to the bounding box
    if city_boundaries is not None:
        gdf_nodes_bb = gdf_nodes[gdf_nodes.intersects(city_boundaries)]
        gdf_edges_bb = gdf_edges[gdf_edges['user_id1'].isin(gdf_nodes_bb['user_id']) & 
                              gdf_edges['user_id2'].isin(gdf_nodes_bb['user_id'])]
        print(f'Filtering to city_boundaries -- nodes: {len(gdf_nodes)}->{len(gdf_nodes_bb)} | edges: {len(gdf_edges)}->{len(gdf_edges_bb)}', flush=True)
        gdf_nodes = gdf_nodes_bb
        gdf_edges = gdf_edges_bb
    return gdf_nodes, gdf_edges

def calculate_social_nodes_density(gdf_nodes, radius_m=400, q=10):
    """
    For every user, estimate the density of their home location by counting other users in a radius
    @param gdf_nodes the dataframe of social network nodes
    @param radius_m: the radius around each node (in meters) to estimate their area density
    @param q: the number of quantiles (10 = deciles)
    @return gdf_nodes: the dataframe with an additional column for local density and one for the quantile 
                       of the density in the empirical density distribution
    """ 
    def _count_spatial_neighbors(row):
        home_location = row['geometry']
        lon = home_location.x
        lat = home_location.y
        buffer = geodesic_point_buffer(lon,lat,radius_m)
        neighbors = len(gdf_nodes[gdf_nodes['geometry'].within(buffer)])
        return neighbors
    gdf_nodes['num_neighbors'] = gdf_nodes.apply(_count_spatial_neighbors, axis=1)
    gdf_nodes['num_neighbors_decile'] = pd.qcut(x=gdf_nodes['num_neighbors'], q=q, duplicates='drop', labels=False)
    return gdf_nodes

def calculate_closest_street_nodes_to_social_nodes(gdf_nodes, G):
    """
    Augments the dataframe with the nodes on the street graph that are closest to the social network nodes
    @param gdf_nodes the dataframe of social network nodes
    @param G: the street graph
    @return gdf_nodes: the dataframe with a new column for the closest street nodes for the home location
    """
    longitudes = gdf_nodes['lon_home']
    latitudes = gdf_nodes['lat_home']
    home_street_nodes = ox.nearest_nodes(G, longitudes, latitudes)
    gdf_nodes['home_street_node'] = home_street_nodes
    return gdf_nodes

def get_shortest_path_length(G_igraph, start_node, dest_node):
    """
    @param G_igraph: the street graph (in iGraph format)
    @param start_node: starting node on the street network
    @param dest_node: destination node on the street network
    @return shortest_path_length: shortest path length in meters
    """
    node_start_index = G_igraph.vs.find(name=start_node).index
    node_dest_index = G_igraph.vs.find(name=dest_node).index
    try:
        shortest_path_length = int(G_igraph.shortest_paths(node_start_index, node_dest_index, weights='length')[0][0])
        return shortest_path_length
    except:
        return None

def get_shortest_path_line(G_igraph, G_gdf_nodes, start_node, dest_node):
    """
    @param G_igraph: the street graph (in iGraph format)
    @param G_gdf_nodes: the GeoPandas dataframe of nodes of the street graph
    @param start_node: starting node on the street network
    @param dest_node: destination node on the street network
    @return route_line: the LineString of the shortest path
    """
    node_start_index = G_igraph.vs.find(name=start_node).index
    node_dest_index = G_igraph.vs.find(name=dest_node).index
    # calculate route with igraph
    route = G_igraph.get_shortest_paths(node_start_index, node_dest_index, weights='length')[0] 
    if route is None:
        return None
    else:
        try:
            # convert iGraph node ids on the route to the original node names
            route = [G_igraph.vs[i]['name'] for i in route]
            route_nodes = G_gdf_nodes.loc[route]
            route_line = LineString(route_nodes['geometry'].tolist())
            return route_line
        except:
            return None

def calculate_social_edges_distances(G_igraph, G_gdf_nodes, gdf_nodes, gdf_edges, ignore_walking=False):
    """
    Augment the edges dataframe with info of density and street network closest nodes
    @param G_igraph: the street graph (in igraph format)
    @param G_gdf_nodes: the nodes of the street network
    @param gdf_nodes: all social network nodes
    @param gdf_edges: all social network edges
    @param ignore_walking: to avoid calculation of walking distance (expensive) - puts walking distance to 0 and geometry to beeline geometry

    """ 
    gdf_edges = pd.merge(gdf_edges, gdf_nodes[['user_id','num_neighbors','num_neighbors_decile','home_street_node']], 
                         left_on='user_id1', right_on='user_id', how='inner')
    del(gdf_edges['user_id'])
    gdf_edges = gdf_edges.rename(columns={'num_neighbors': 'num_neighbors1', 
                                          'num_neighbors_decile': 'num_neighbors_decile1',
                                          'home_street_node': 'home_street_node1'})
    gdf_edges = pd.merge(gdf_edges, gdf_nodes[['user_id','num_neighbors','num_neighbors_decile','home_street_node']], 
                         left_on='user_id2', right_on='user_id', how='inner')
    del(gdf_edges['user_id'])
    gdf_edges = gdf_edges.rename(columns={'num_neighbors': 'num_neighbors2', 
                                          'num_neighbors_decile': 'num_neighbors_decile2',
                                          'home_street_node': 'home_street_node2'})
    
    # calculate the beeline distance of edges
    gdf_edges['distance_beeline_m'] = gdf_edges.apply(lambda x : int(great_circle((x['lat_home1'], x['lon_home1']), (x['lat_home2'], x['lon_home2'])).m), axis=1)
    # zero-distance is not meaningful remove it
    print(f'removing {len(gdf_edges[gdf_edges.distance_m == 0])} edges with zero-distance')
    gdf_edges = gdf_edges[gdf_edges['distance_beeline_m'] > 0]
    if not ignore_walking:
        # calculate the geometry of the shortes path
        gdf_edges['shortest_path_geometry'] = gdf_edges.apply(lambda x : get_shortest_path_line(G_igraph, G_gdf_nodes, x['home_street_node1'], x['home_street_node2']), axis=1)
        # calculate the walking distance of edges
        gdf_edges['distance_walking_m'] = gdf_edges.apply(lambda x : get_shortest_path_length(G_igraph, x['home_street_node1'], x['home_street_node2']), axis=1)
    else:
        gdf_edges['shortest_path_geometry'] = gdf_edges['geometry']
        gdf_edges['distance_walking_m'] = 0

    return gdf_edges

def keep_reciprocal(df, field1, field2):
    """
    Returns a dataframe containing only reciprocal entries in field1 and field2 
    (i.e., all those x,y entries for which there exists a y,x pair in the dataframe)
    @param df: a dataframe
    @param field1: the first field
    @param field2: the second field
    @return df_reciprocal: dataframe containing only reciprocal entries
    """
    df_reciprocal_pairs = df[[field1, field2]]
    df_reciprocal_pairs['sorted1'] = df_reciprocal_pairs.apply(lambda x: min(x[field1], x[field2]), axis=1)
    df_reciprocal_pairs['sorted2'] = df_reciprocal_pairs.apply(lambda x: max(x[field1], x[field2]), axis=1)
    df_reciprocal_pairs.sort_values(by=['sorted1', 'sorted2'])
    df_reciprocal_pairs = df_reciprocal_pairs.groupby(['sorted1', 'sorted2']).count().reset_index()
    df_reciprocal_pairs = df_reciprocal_pairs[df_reciprocal_pairs[field1] == 2]
    df_reciprocal_pairs = df_reciprocal_pairs[['sorted1', 'sorted2']]
    df_reciprocal_pairs.columns = [field1, field2]
    df2 = pd.DataFrame({field1:df_reciprocal_pairs[field2], field2:df_reciprocal_pairs[field1]})
    df_reciprocal_pairs = pd.concat([df_reciprocal_pairs, df2])
    return pd.merge(df, df_reciprocal_pairs)

def get_filename_CM_edges(
    cbsacode,
    distance_type,
    filter_by_density,
    random_seed,
    output_folder = "../data/social_networks_null/"
):
    '''helper function for configuration_model_spatial_constraint() filename lengthiness'''
    return output_folder + f'{cbsacode}_gdf_edges_null_{distance_type}_densityfilter-{filter_by_density}_s{random_seed}.pq'

def get_filename_CM_logs(
    cbsacode,
    distance_type,
    filter_by_density,
    random_seed,
    output_folder = "../data/error_logs/"
):
    '''helper function for configuration_model_spatial_constraint() filename lengthiness'''
    return output_folder + f'{cbsacode}_error_logs_{distance_type}_densityfilter-{filter_by_density}_s{random_seed}.pickle'

def configuration_model_spatial_constraint(
        gdf_edges, 
        gdf_nodes, 
        G, 
        G_igraph, 
        G_gdf_nodes, 
        cbsacode, 
        city_boundaries,
        distance_type = 'beeline', 
        geo_buffer = 150, 
        geo_buffer_step = 50, 
        geo_buffer_step_multiplier = 1.03, 
        filter_by_density=False, 
        calculate_walking_on_selected=False,
        test_mode=False, 
        undirected=False, 
        random_seed=42,
        outfolder_nullmodel="../data/social_networks_null/",
        outfolder_errorlogs="../data/error_logs/",
        outfolder_plots="../plots/"
        ):
    
    """
    Build a randomly rewired network that preserves degree sequence (configuration model),
    distances, and preferential attachment for density
    
    @param gdf_edges: dataframe with social edges
    @param gdf_nodes: dataframe with nodes
    @param G: the street graph
    @param G_igraph: the street graph (in igraph format)
    @param cbsacode: cbsa code of the region of interest
    @param distance_type: the type of distance to be calculated between nodes when selecting them: ['beeline' | 'walking']
    @param geo_buffer: initial geographic buffer 
    @param geo_buffer_step: how many meters is added to the buffer if no candidate is found
    @param geo_buffer_step_multiplier: multiplier that slowly increases the width of the geo buffer step
    @param buffer_initial: initial fraction of buffer enlargement
    @param filter_by_density: toggles filter by spatial density
    @param calculate_walking_on_selected: calculate walking distance and path on the selected node
    @param test_mode: toggle for testing (just one iteration, print a map of a rewired edge)
    @param undirected: to handle undirected networks (patent dataset)
    @param random_seed: the random seed for the randomized model
    """

    # select desired type of distance 
    _supported_distance_types = ['beeline', 'walking']
    if distance_type == 'beeline':
        distance_column = 'distance_beeline_m'
        edge_geometry_column = 'geometry'
    elif distance_type == 'walking':
        distance_column = 'distance_walking_m'
        edge_geometry_column = 'shortest_path_geometry'
        # recalculates walking distances from path
        geod = Geod(ellps="WGS84")
        gdf_edges[distance_column] = gdf_edges[edge_geometry_column].apply(lambda x : round(geod.geometry_length(x)))
        print('distances from path recalculated', flush=True)
    else:
        raise Exception(f'Supported distance types are {_supported_distance_types}')
    #gdf_edges[f'distance_{distance_type}_m_approx'] = gdf_edges[distance_column].apply(lambda x: ((x//100)*100)+100)
    
    res = [] # result will be saved here
    dist_error_log = [] # log the absolute errors in distance
    dist_error_perc_log = [] # log the % errors in distance
    tolerance_log = [] # log the absolute errors in density deciles
    cnt = 0 # iteration count
    cache_hits = 0 # walking distance polygon cache
    
    isochrones_cache = {} # chaches the isochrone buffers for nodes
    
    # important step of re-indexing, it is crucial for the configuration model
    #gdf_nodes = gdf_nodes.reset_index()
    gdf_nodes =  gdf_nodes.set_index('user_id')
    gdf_edges = gdf_edges.reset_index()    

    # get list of stubs (incoming links) for the configuration model
    df_in_stubs = gdf_edges[['user_id2', 'cbsacode']].groupby('user_id2').count().reset_index()
    df_in_stubs.columns = ['user_id', 'in_stubs']
    if undirected: #adds also out degree
        df_out_stubs = gdf_edges[['user_id1', 'cbsacode']].groupby('user_id1').count().reset_index()
        df_out_stubs.columns = ['user_id', 'in_stubs']
        df_in_stubs = pd.concat([df_in_stubs, df_out_stubs]).groupby('user_id').sum().reset_index()
    df_in_stubs = df_in_stubs.set_index('user_id')

    print('==', len(df_in_stubs), len(gdf_edges), flush=True)

    # restrict to short edges when calculating walking
    if distance_type == 'walking':
        gdf_edges = gdf_edges[gdf_edges['distance_beeline_m'] <= 12000] 

    for idx, row in tqdm(gdf_edges.sample(frac=1, random_state=random_seed).iterrows()):
        if distance_type == 'walking':
            if cnt >= 10000:
                break
        cnt+=1
        if cnt%1000 == 0:
            print(f'{cnt} edges parsed ({cache_hits} cache hits)', flush=True)
            print('---', len(df_in_stubs), df_in_stubs['in_stubs'].sum(), flush=True) ############*********
        # get information about current nodes and edge
        user_id1 = row['user_id1']
        user_id2_original = row['user_id2']
        lon_home1 = row['lon_home1'] 
        lat_home1 = row['lat_home1']
        lon_home2_original = row['lon_home2'] 
        lat_home2_original = row['lat_home2']
        location_node_origin = Point(lon_home1, lat_home1)
        original_edge_geom = row[edge_geometry_column]
        dist_beeline = row[f'distance_beeline_m']
        user_id1_G_node = row['home_street_node1']
        user_id2_original_G_node = row['home_street_node2']
        dist_walking = row[f'distance_walking_m']    
        if distance_type == 'beeline':
            dist_original = dist_beeline
        if distance_type == 'walking':
            dist_original = dist_walking
        density = row['num_neighbors_decile2']
        #print(f'orginalrow = {user_id1}, {row["user_id2"]} ({lon_home1}, {lat_home1}) -> ({lon_home2}, {lat_home2}), {dist_exact}, {original_edge_geom}')  

        # consider only remaining nodes with available stubs
        df_available_nodes = gdf_nodes[gdf_nodes.index.isin(df_in_stubs.index)]
        #print('---', len(df_available_nodes), flush=True) ############*********

        # filter by distance
        df_candidates = []
        geo_buffer = 0
        if distance_type == 'walking':
            geo_buffer += geo_buffer_step
        attempts = 10
        j = 0 ############*********
        while len(df_candidates) == 0 and j<attempts: ############*********
            j +=1 ############*********
            if distance_type == 'beeline':
                try:
                    # draw a torus around the home location, with diameter similar to the distance to actual friend
                    buffer_outer = geodesic_point_buffer(lon_home1, lat_home1, dist_original + geo_buffer)
                    buffer_inner = geodesic_point_buffer(lon_home1, lat_home1, max(10, dist_original - geo_buffer))
                    #buffer_ring = Polygon(buffer_outer.exterior.coords, [buffer_inner.exterior.coords])                
                    #if geo_buffer >= 50:
                    #    print('buf', geo_buffer)
                except: # this should happen very rarely
                    buffer_outer = None 
                    buffer_inner = None
            elif distance_type == 'walking':
                try:
                    buffer_outer, buffer_inner = isochrones_cache[(home_street_node1, (dist_original//100)*100)]
                    cache_hits+=1
                except:
                    outer_subgraph = nx.ego_graph(G, user_id1_G_node, radius=dist_original+geo_buffer, distance='length')
                    outer_node_points = [Point((data['x'], data['y'])) for node, data in outer_subgraph.nodes(data=True)]
                    buffer_outer = gpd.GeoSeries(outer_node_points).unary_union.convex_hull
                    inner_subgraph = nx.ego_graph(G, user_id1_G_node, radius=max(10, dist_original-geo_buffer), distance='length')
                    inner_node_points = [Point((data['x'], data['y'])) for node, data in inner_subgraph.nodes(data=True)]
                    buffer_inner = gpd.GeoSeries(inner_node_points).unary_union.convex_hull
                    #buffer_ring = Polygon(buffer_outer.exterior.coords, [buffer_inner.exterior.coords])                
                    isochrones_cache[(user_id1_G_node,dist_original)] = [buffer_outer, buffer_inner]
                    #print(len(outer_subgraph), len(inner_subgraph), user_id1_G_node, dist_original)
            else:
                raise Exception(f'Supported distance types are {_supported_distance_types}')
            # get candidates inside the torus
            if buffer_outer is not None:
                df_candidates = df_available_nodes[df_available_nodes['geometry'].within(buffer_outer) 
                                                   & ~df_available_nodes['geometry'].within(buffer_inner)]
            else:
                print(f'{cbsacode} >>> unable to get a buffer, selecting random candidate')
                df_candidates = df_available_nodes
            geo_buffer += geo_buffer_step
            geo_buffer *= geo_buffer_step_multiplier
            if distance_type == 'walking':
                geo_buffer *= geo_buffer_step_multiplier # do it twice for faster convergence
        
        if len(df_candidates) == 0 : #after all attempts nothing could be found
            print(f'**, {geo_buffer}', flush=True) ############*********
            continue

        # filter by density OLD
        # if filter_by_density:
        #     tolerance = tolerance_initial
        #     df_candidates_density = []
        #     while len(df_candidates_density) == 0:
        #         df_candidates_density = df_candidates[abs(df_candidates['num_neighbors_decile'] - density) <= tolerance]
        #         tolerance += tolerance_step
        #         #print('>>>', len(df_candidates_density))
        #         #if tolerance >= 2:
        #         #    print('tol', tolerance)
        #     tolerance_log.append(tolerance-1)
        # else:
        #     df_candidates_density = df_candidates

        # select a random candidate among the candidates remaining after the filters
        if filter_by_density: # selection proportional to local density of the destination node
            df_candidates['cum'] = df_candidates['num_neighbors'].cumsum()
            min_v = df_candidates['cum'].min()
            max_v = df_candidates['cum'].max()
            v = random.randint(min_v, max_v)
            selected_user = pd.merge(df_in_stubs, df_candidates[df_candidates['cum']>=v].head(1), on='user_id')
        else: # selection completely at random
            selected_user = pd.merge(df_in_stubs, df_candidates.sample(1, random_state=random_seed), on='user_id')

        # update the available stubs
        ######selected_user_id = selected_user['user_id']
        selected_user_id = selected_user.index
        # df_in_stubs.to_csv('xxxx',index=True)
        # print('------', selected_user_id, user_id1, flush=True)
        # print(1, df_in_stubs.loc[selected_user_id], flush=True)
        # print(2, df_in_stubs.loc[selected_user_id]['in_stubs'], flush=True)
        # print(3, df_in_stubs.loc[selected_user_id]['in_stubs'].values, flush=True)
        # print(4, df_in_stubs.loc[selected_user_id]['in_stubs'].values[0], flush=True)
        num_stubs = df_in_stubs.loc[selected_user_id]['in_stubs'].values[0]
        if num_stubs == 1:
            df_in_stubs = df_in_stubs.drop(selected_user_id)
        else:
            df_in_stubs.loc[selected_user_id, 'in_stubs'] = num_stubs-1
            #df_in_stubs.loc[selected_user_id]['in_stubs'] = num_stubs-1
        # in the undirected case, the source node also loses its stub
        if undirected:
            # print(1, df_in_stubs.loc[user_id1], flush=True)
            # print(2, df_in_stubs.loc[user_id1]['in_stubs'], flush=True)
            # print(3, df_in_stubs.loc[user_id1]['in_stubs'].values, flush=True)
            # print(4, df_in_stubs.loc[user_id1]['in_stubs'].values[0], flush=True)
            num_stubs = df_in_stubs.loc[user_id1]['in_stubs']
            if num_stubs == 1:
                df_in_stubs = df_in_stubs.drop(user_id1)
            else:
                df_in_stubs.loc[user_id1, 'in_stubs'] = num_stubs-1

        #print(f'remaining stubs: {df_in_stubs.in_stubs.sum()}')

        # get info on new edge and calculate exact new distance 
        user_id2 = selected_user_id._data[0]
        lat_home2, lon_home2, location_node_destination = selected_user[['lat_home','lon_home','geometry']].values[0]
        
        newdist_beeline = int(distance((lat_home1,lon_home1),(lat_home2,lon_home2)).m)
        user_id2_G_node = None
        newdist_walking = None
        geometry_walking = None 
        
        geometry = LineString([location_node_origin, location_node_destination])

        # calculate new walking distance and path
        if calculate_walking_on_selected or distance_type == 'walking':
            try:
                geometry_walking = get_shortest_path_line(G_igraph, G_gdf_nodes, user_id1_G_node, user_id2_G_node)
                user_id2_G_node = gdf_nodes.loc[user_id2]['home_street_node']
                newdist_walking = get_shortest_path_length(G_igraph, user_id1_G_node, user_id2_G_node)
            except:
                geometry_walking = None 
                newdist_walking = None
        
        if distance_type == 'beeline':
            dist_err = abs(dist_beeline - newdist_beeline)
            dist_err_perc = dist_err/dist_beeline
        elif distance_type == 'walking':
            if newdist_walking:
                dist_err = abs(dist_walking - newdist_walking)
            else:
                dist_err = 0 # this should not happen
            dist_err_perc = dist_err/dist_walking
        dist_error_log.append(dist_err)
        dist_error_perc_log.append(dist_err_perc)
        
        # add new edge to the result list
        res.append([cbsacode, user_id1, user_id2, user_id2_original, lat_home1, lon_home1, 
                    lat_home2, lon_home2, lat_home2_original, lon_home2_original,
                    user_id1_G_node, user_id2_G_node, user_id2_original_G_node,
                    dist_beeline, dist_walking, newdist_beeline, newdist_walking, geometry, geometry_walking])
        
        # enable just for testing
        if test_mode:
            gdf_edges_null = gpd.GeoDataFrame(res, columns=['cbsacode', 'user_id1', 'user_id2', 'user_id2_original', 'lat_home1','lon_home1', 
                                                            'lat_home2', 'lon_home2', 'lat_home2_original', 'lon_home2_original', 'home_street_node1',
                                                            'home_street_node2', 'home_street_node2_original', 'distance_beeline_m_orginal', 
                                                            'distance_walking_m_original', 'distance_beeline_m', 'distance_walking_m', 'geometry', 'geometry_walking'],
                                     geometry="geometry", crs="epsg:4326")
            fig, ax = plt.subplots(1,1, figsize=(10,10))
            ax.plot(*city_boundaries.exterior.xy)
            #gdf_edges_city.plot(lw=0.1, color='orange', alpha=0.5, ax=ax)
            gdf_nodes.plot(markersize=2, color='yellow', alpha=0.3, ax=ax)
            ax.plot(*original_edge_geom.xy, lw=3, color='white')
            ax.scatter([lon_home1],[lat_home1], s=60, c='white')
            ax.scatter([lon_home2],[lat_home2], s=60, c='red')
            gdf_edges_null.plot(lw=3, color='red', ax=ax)
            ax.plot(*buffer_outer.exterior.xy, c='white')
            ax.plot(*buffer_inner.exterior.xy, c='white')
            #gdf_edges.plot(lw=0.1, color='white', alpha=0.1, ax=ax)
            ax.axis('off')
            plt.tight_layout()
            fig.patch.set_facecolor('black')
            plt.savefig(outfolder_plots + f'{cbsacode}_map_test.png', dpi=600)
    #         candidates = df_candidates['geometry'].values
    #         fig, ax = plt.subplots(1,1, figsize=(5,5))
    #         plt.scatter([lon_home1], [lat_home1], c='purple')
    #         plt.plot(*buffer_outer.exterior.xy, c='red', linewidth=3)
    #         plt.plot(*buffer_inner.exterior.xy, c='red', linewidth=3)        
    #         xs = [point.x for point in candidates]
    #         ys = [point.y for point in candidates]
    #         plt.scatter(xs, ys)
            if cnt>0:
                return None
        # if cnt > 1000:
        #     break
            
    # create a dataframe with the randomized edges
    gdf_edges_null = gpd.GeoDataFrame(res, columns=['cbsacode', 'user_id1', 'user_id2', 'user_id2_original', 'lat_home1','lon_home1', 
                                                    'lat_home2', 'lon_home2', 'lat_home2_original', 'lon_home2_original', 'home_street_node1',
                                                    'home_street_node2', 'home_street_node2_original', 'distance_beeline_m_orginal', 
                                                    'distance_walking_m_original', 'distance_beeline_m', 'distance_walking_m', 'geometry', 'geometry_walking'],
                             geometry="geometry", crs="epsg:4326")
    # join with original edges to get original geometry
    df_join = gdf_edges[['user_id1', 'user_id2', edge_geometry_column]]
    df_join.columns = ['user_id1', 'user_id2_original', 'geometry_original']
    gdf_edges_null = pd.merge(gdf_edges_null, df_join, on = ['user_id1', 'user_id2_original'])
    gdf_edges_null['geometry_walking'] = gpd.geoseries.GeoSeries(gdf_edges_null['geometry_walking'])
    
    # saving results on file
    outfile_name = get_filename_CM_edges(cbsacode, distance_type, filter_by_density, random_seed, outfolder_nullmodel)
    gdf_edges_null.to_parquet(outfile_name)
    
    errs = {'dist_error_log': dist_error_log, 
        'dist_error_perc_log': dist_error_perc_log,
        'tolerance_log':tolerance_log}

    name_fout = get_filename_CM_logs(cbsacode, distance_type, filter_by_density, random_seed, outfolder_errorlogs) 
    with open(name_fout, 'wb') as handle:
            pickle.dump(errs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    # printing a couple of check measures
    print(f'remaining stubs: {df_in_stubs.in_stubs.sum()}')
    print(f'cache hits: {cache_hits}')
    return gdf_edges_null, dist_error_log, dist_error_perc_log, tolerance_log

def sanity_degree_sequence(
        gdf_edges, 
        gdf_edges_null, 
        cbsacode, 
        plot=True, 
        save=True, 
        outfolder_plots = "../plots/"
        ):
    """
    Check that the degree sequence of the null model is identical to real data
    @param gdf_edges: dataframe with social edges
    @param gdf_edges_null: dataframe with social edges (reshuffled)
    @param cbsacode: cbsa code of the region of interest
    """ 
    res = {'indegree': None, 'outdegree': None}
    res['outdegree'] = ks_2samp(gdf_edges.groupby('user_id1').count()['cbsacode'], 
                               gdf_edges_null.groupby('user_id1').count()['cbsacode'])
    res['indegree'] = ks_2samp(gdf_edges.groupby('user_id2').count()['cbsacode'], 
                               gdf_edges_null.groupby('user_id2').count()['cbsacode'])
    if plot:
        fig, axs = plt.subplots(1,2, figsize=(15,5))
        ax1, ax2 = axs
        gdf_edges.groupby('user_id1').count()['cbsacode'].hist(bins=100, ax=ax1, alpha=0.5)
        gdf_edges_null.groupby('user_id1').count()['cbsacode'].hist(bins=100, ax=ax1, alpha=0.5)
        ax1.set_xlabel('outdegree')
        ax1.set_ylabel('frequency')

        gdf_edges.groupby('user_id2').count()['cbsacode'].hist(bins=100, ax=ax2, alpha=0.5)
        gdf_edges_null.groupby('user_id2').count()['cbsacode'].hist(bins=100, ax=ax2, alpha=0.5)
        ax2.set_xlabel('indegree')
        plt.tight_layout()
        if save:
            plt.savefig(outfolder_plots + f'{cbsacode}_sanity1.png', dpi=600)
    return res

def sanity_errors(
        dist_error_log, 
        dist_error_perc_log, 
        tolerance_log, 
        cbsacode, 
        plot=True, 
        save=True, 
        outfolder_plots = "../plots/"
        ):
    """
    Measure errors in terms of distances and density compared to the original data
    @param dist_error_log: list of distance errors (meters)
    @param dist_error_perc_log: list of distance errors (%)
    @param tolerance_log: list of distance errors (%)
    @param cbsacode: cbsa code of the region of interest
    """ 
    fig, axs = plt.subplots(1,3, figsize=(15,4))
    ax1, ax2, ax3 = axs
    
    res = {'error_meters': (np.mean(dist_error_log), np.median(dist_error_log) ),
           'error_perc': (np.mean(dist_error_perc_log), np.median(dist_error_perc_log) ),
           'error_tolerance': (np.mean(tolerance_log), np.median(tolerance_log) ),
          }
    
    if plot:
        y,x = np.histogram(dist_error_log, bins = 2000)
        ax1.scatter(x[:-1], y, c='#E84893')
        ax1.set_xlabel('error on distance (m)')
        ax1.set_ylabel('# edges')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax1.set_xlim(0,10000)

        y,x = np.histogram(dist_error_perc_log, bins = 3000) 
        ax2.scatter(x[:-1], y, c='#E84893')
        ax2.set_xlabel('% error on distance')
        ax2.set_ylabel('# edges')
        #ax2.set_xlim(0,100)
        ax2.set_xscale('log')
        ax2.set_yscale('log')

        y,x =np.histogram(tolerance_log, bins = 10)
        ax3.scatter(x[:-1], y, c='#E84893')
        ax3.set_xlabel('tolerance on target node density')
        ax3.set_ylabel('# edges')
        plt.tight_layout()
        if save:
            plt.savefig(outfolder_plots + f'{cbsacode}_sanity2.png', dpi=600)
    return res

def sanity_distance_distribution(
        gdf_edges_null, 
        cbsacode, 
        plot=True, 
        save=True,
        outfolder_plots = "../plots/"
        ):
    """
    Measure difference in the distance distributions
    @param gdf_edges: dataframe with social edges
    @param gdf_edges_null: dataframe with social edges (reshuffled)
    @param cbsacode: cbsa code of the region of interest
    """ 
    res = {}
    for dist_type in ['beeline', 'walking']:
        res[dist_type] = ks_2samp(gdf_edges_null[f'distance_{dist_type}_m_orginal'], gdf_edges_null[f'distance_{dist_type}_m'])
        if plot:
            fig, ax = plt.subplots(1,1, figsize=(15,5))
            gdf_edges_null[f'distance_{dist_type}_m_orginal'].hist(bins=200, ax=ax, alpha=0.5, label='original')
            gdf_edges_null[f'distance_{dist_type}_m'].hist(bins=200, ax=ax, alpha=0.5, label='null')
            ax.set_ylabel('number of edges')
            ax.set_xlabel('distance (m)')
            plt.legend()
            plt.tight_layout()
            if save:
                plt.savefig(outfolder_plots + f'{cbsacode}_sanity3.png', dpi=600)
        return res

def plot_street_network(
        gdf_edges_city, 
        cbsacode, 
        outfolder_plots = "../plots/"
    ):
    
    # plotting the street network
    street_types = ['motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'residential']
    colors = ['#C33478', '#C33478', '#CF6F3B', '#90b7f2', '#90b7f2', 'gray']
    width = [3, 3, 2, 0.8, 0.8, 0.3]
    alpha = [0.5, 0.5, 0.5, 0.5, 0.5, 0.2]
    i = len(street_types)
    fig, ax = plt.subplots(1,1, figsize=(15,15))
    for street_type, c , w, a  in zip(street_types, colors, width, alpha):
        df_selection = gdf_edges_city[gdf_edges_city['highway']==street_type]
        if len(df_selection)>0: # plot street types if they exist 
            df_selection = gpd.GeoDataFrame(df_selection, geometry=df_selection['geometry'], crs='epsg:4326')
            df_selection.plot(lw=w, color=c, alpha=a, ax=ax, label=street_type, zorder = i)
            i -= 1
    fig.patch.set_facecolor('#f3ebd5')
    ax.legend()
    #gdf_nodes_bb.plot(markersize=2, color='yellow', alpha=0.3, ax=ax)
    #gdf_edges.plot(lw=0.1, color='white', alpha=0.1, ax=ax)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(outfolder_plots + f'{cbsacode}_streetmap.png', dpi=800)

def intersect_social_and_street(gdf_edges, gdf_edges_null, gdf_edges_city, distance_type='beeline'):

    # to cast types
    dtypes = {'cbsacode': str, 'user_id1': int, 'user_id2': int, 'lat_home1': float, 'lon_home1': float,
                'lat_home2': float, 'lon_home2': float, 'distance_beeline_m': int, 'distance_walking_m': int, 
                'osmid': str, 'highway': str, 'maxspeed': str, 'oneway': bool, 'reversed': bool,
                'length': float, 'lanes': str, 'service': str, 'bridge': str, 'name': str, 'width': str, 
                'ref': str, 'access': str, 'tunnel': str, 'junction': str}


    # # DASK of social edges
    # ddf_edges = gdf_edges.copy()
    # ddf_edges = ddf_edges.reset_index()
    # ddf_edges = dask_geopandas.from_geopandas(ddf_edges, npartitions=16)

    # # DASK of null model social edges
    # ddf_edges_null = gdf_edges_null.copy()
    # ddf_edges_null = ddf_edges_null.reset_index()
    # ddf_edges_null = dask_geopandas.from_geopandas(gdf_edges_null, npartitions=16)
    # if distance_type == 'walking':
    #     ddf_edges_null = ddf_edges_null.dropna()

    # # DASK of street edges
    # ddf_edges_city = gdf_edges_city.copy()
    # ddf_edges_city = ddf_edges_city.reset_index()
    # ddf_edges_city = dask_geopandas.from_geopandas(ddf_edges_city, npartitions=16)

    # try:
    #     print('intersecting with dask', flush=True)
    # fields = ['cbsacode', 'user_id1', 'user_id2', 'lat_home1', 'lon_home1', 'lat_home2', 'lon_home2', 
    #             'home_street_node1', 'home_street_node2', 'distance_beeline_m', 'distance_walking_m', 
    #             'index_right', 'u', 'v', 'key', 'osmid', 'highway', 'maxspeed', 'oneway', 'reversed', 
    #             'length', 'lanes', 'service', 'bridge', 'name', 'width', 'ref', 'access', 'tunnel', 
    #             'junction']
    #     ddf = ddf_edges.sjoin(ddf_edges_city, how='inner')
    #     print('dask succeeded!')
    #     ddf = ddf[fields]
    #     ddf_null = ddf_edges_null.sjoin(ddf_edges_city, how='inner')
    #     ddf_null = ddf_null[fields]
    #     for attr in dtypes:
    #         ddf[attr] = ddf[attr].astype(dtypes[attr])
    #         ddf_null[attr] = ddf_null[attr].astype(dtypes[attr])

    #     ddf = ddf.compute()
    #     ddf_null = ddf_null.compute()
        
    #     # reconvert to geopandas GeoDataFrame
    #     gdf = dask_geopandas.from_dask_dataframe(ddf)
    #     gdf_null = dask_geopandas.from_dask_dataframe(ddf_null)
   
    # except:
    #     print("dask failed.")
    #     print("intersecting with geopandas")
    gdf = gpd.sjoin(
        gdf_edges, 
        gdf_edges_city, 
        lsuffix='social', 
        rsuffix='street', 
        how='inner')
    gdf_null = gpd.sjoin(
        gdf_edges_null, 
        gdf_edges_city, 
        lsuffix='social', 
        rsuffix='street', 
        how='inner'
        )
    for attr in dtypes:
        if attr in gdf.columns:
            gdf[attr] = gdf[attr].astype(dtypes[attr])
            gdf_null[attr] = gdf_null[attr].astype(dtypes[attr])
    
    return gdf, gdf_null

def get_random_street_network(G, G_igraph, gdf_nodes_city, street_lengths):
    """
    @param G: the street graph
    @param G_igraph: the street graph (in igraph format)
    @param gdf_nodes_city: the nodes of the street graph
    @param street_lengths: a list of street lengths to be reproduced
    @param min_residual_length: if the sum of lengths of remaining streets is below this threshold, 
                                only one street with length ~equal to the remaining sum is created
    @param tolerance: street lengths are matched with some % tolerance (value in [0,1])
    @return random_streets, random_street_lengths: lists of street LineStrings and their lengths
    """
    street_lengths = copy.deepcopy(street_lengths)
    tot_length = sum(street_lengths)
    print(f'#streets={len(street_lengths)}, total_length={tot_length}', flush=True)
    random_streets = []
    random_street_lengths = []
    while tot_length > 500: # if remaining street length is negligible, stop
        # generate a random pair of street nodes
        start_node_id = dest_node_id = 0
        while start_node_id == dest_node_id:
            start_node_id = choice(list(G.nodes))
            dest_node_id = choice(list(G.nodes))
        start_node = G.nodes[start_node_id]
        dest_node = G.nodes[dest_node_id]
        # calculate distance between those nodes
        d = int(great_circle((start_node['y'],start_node['x']), (dest_node['y'], dest_node['x'])).m)
        if d > max(street_lengths) or d > tot_length:
            continue
        # trace a path between the two points
        random_path_length = get_shortest_path_length(G_igraph, start_node_id, dest_node_id)
        if random_path_length is None or random_path_length <= 0:
            continue
        tot_length = tot_length - random_path_length
        random_path_line = get_shortest_path_line(G_igraph, gdf_nodes_city, start_node_id, dest_node_id)
        random_streets.append(random_path_line)
        random_street_lengths.append(random_path_length)
        #print(f'tot_remain={tot_length}, tot_rand={sum(random_street_lengths)} | rand_bee={d}, rand_walk={random_path_length}')            
    return random_streets, random_street_lengths