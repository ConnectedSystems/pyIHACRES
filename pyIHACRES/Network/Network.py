import yaml

from .DamNode import DamNode
from .StreamNode import StreamNode


class Network(object):

    """Represent a stream network"""

    def __init__(self, network_details, first_node):
        self.node_types = {
            'StreamNode': StreamNode,
            'DamNode': DamNode
        }

        network = {}
        for node_id in network_details:
            self.construct_nodes(node_id, network_details, network)
        # End for

        self.network_details = network_details
        self.network = network
        self.first_node = first_node
        self.ts = 0
    # End init()

    @staticmethod
    def load_network(fn, first_node='node_1'):
        with open(fn, 'r') as network_config:
            network_details = yaml.load(network_config, Loader=yaml.FullLoader)

        return Network(network_details, first_node)
    # End load_network()

    def construct_nodes(self, node_id, details, network):
        """Recursively construct nodes.

        :param node_id: str, node identifier
        :param details: dict, of all node parameters
        :param network: dict, resulting stream network

        :returns: object, StreamNode or DamNode, or None
        """
        if node_id in network:
            return network[node_id]

        node_details = details.get(node_id, None)
        if node_details:
            node_details.update({'node_id': node_id})
        else:
            return None
        # End if

        params = {k: v for k, v in node_details.items() if k != 'node_type'}
        this_node = self.node_types[node_details['node_type']](**params)
        network[node_id] = this_node

        this_node.next_node = self.construct_nodes(this_node.next_node, details, network)

        # Construct previous nodes, if any
        for prev_id in this_node.prev_node:
            this_node.prev_node[prev_id] = self.construct_nodes(prev_id, details, network)
        # End for

        return this_node

    # End construct_nodes()

    def run_timestep(self, rain, et, irrig_ext, extractions):
        """Run model for a time step.

        Each parameter should be a dict for each node in the network.

        :param rain: dict, rainfall (in mm) for each node
        :param et: dict, evapotranspiration (in mm) for each node
        :param irrig_ext: dict, irrigation extractions (in ML) for each node
        :param extractions: dict, other extractions (in ML) for each node
        """
        node = self.network[self.first_node]
        while node is not None:
            node = self.run_node(node)
        # End while

    # End run_timestep()

    def reset(self):
        tmp = Network(self.network_details, self.first_node)
        self.network = tmp.network
    # End reset()
