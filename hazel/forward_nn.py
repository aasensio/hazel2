import numpy as np
import glob
import torch
import torch_geometric.data
import torch_geometric.loader
import hazel.graphnet as graphnet
from sklearn import neighbors
import logging

class Dataset(torch.utils.data.Dataset):
    def __init__(self, hyperparameters, tau_all, ne_all, vturb_all, T_all, vlos_all):
        """
        Dataset for the depth stratification
        """
        super(Dataset, self).__init__()


        # Now we need to define the graphs for each one of the computed models
        # The graph will connect all points at certain distance. We define this distance
        # as integer indices, so that we make sure that nodes are connected to the neighbors
        self.n_training = len(T_all)

        # Initialize the graph information
        self.edge_index = [None] * self.n_training
        self.nodes = [None] * self.n_training
        self.edges = [None] * self.n_training
        self.u = [None] * self.n_training

        # Loop over all training examples
        for i in range(self.n_training):

            num_nodes = len(tau_all[i])
            index_cmass = np.zeros((num_nodes, 1))

            index_cmass[:, 0] = np.arange(num_nodes)  # self.cmass_all

            # Build the KDTree
            self.tree = neighbors.KDTree(index_cmass)

            # Get neighbors
            receivers_list = self.tree.query_radius(index_cmass, r=1)

            senders = np.repeat(range(num_nodes), [len(a) for a in receivers_list])
            receivers = np.concatenate(receivers_list, axis=0)            

            # Mask self edges
            mask = senders != receivers            

            # Transform senders and receivers to tensors
            senders = torch.tensor(senders[mask].astype('long'))
            receivers = torch.tensor(receivers[mask].astype('long'))
            

            # Define the graph for this model by using the sender/receiver information
            self.edge_index[i] = torch.cat([senders[None, :], receivers[None, :]], dim=0)

            n_edges = self.edge_index[i].shape[1]

            # Now define the nodes. For the moment we use only one quantity, the log10(T)
            node_input_size = hyperparameters['node_input_size']
            self.nodes[i] = np.zeros((num_nodes, node_input_size))
            self.nodes[i][:, 0] = np.log10(T_all[i])
            self.nodes[i][:, 1] = np.log10(tau_all[i])            
            self.nodes[i][:, 2] = np.log10(ne_all[i])
            self.nodes[i][:, 3] = vturb_all[i] / 1e3
            self.nodes[i][:, 4] = vlos_all[i] / 1e3

            # We use two quantities for the information encoded on the edges: log(column mass) and log(tau)
            edge_input_size = hyperparameters['edge_input_size']
            self.edges[i] = np.zeros((n_edges, edge_input_size))
            
            if edge_input_size == 1:
                tau0 = np.log10(tau_all[i][self.edge_index[i][0, :]])
                tau1 = np.log10(tau_all[i][self.edge_index[i][1, :]])
                self.edges[i][:, 0] = (tau0 - tau1)
            else:
                raise ValueError("Incompatible edge input size")

            # We don't use at the moment any global property of the graph, so we set it to zero.
            self.u[i] = np.zeros((1, 1))
            # self.u[i][0, :] = np.array([np.log10(self.eps_all[i][0, 0]), np.log10(self.ratio_all[i][0, 0])], dtype=np.float32)

            # Finally, all information is transformed to float32 tensors
            self.nodes[i] = torch.tensor(self.nodes[i].astype('float32'))
            self.edges[i] = torch.tensor(self.edges[i].astype('float32'))
            self.u[i] = torch.tensor(self.u[i].astype('float32'))
        
    def __getitem__(self, index):

        # When we are asked to return the information of a graph, we encode
        # it in a Data class. Batches in graphs work slightly different than
        # in more classical situations. Since we have the connectivity of each
        # graph, batches are built by generating a big graph containing all
        # graphs of the batch.
        node = self.nodes[index]
        edge_attr = self.edges[index]        
        u = self.u[index]
        edge_index = self.edge_index[index]

        data = torch_geometric.data.Data(x=node, edge_index=edge_index, edge_attr=edge_attr, y=None, u=u)

        return data

    def __len__(self):
        return self.n_training

    def __call__(self, index):
        return self.cmass_all[index], self.tau_all[index], self.vturb_all[index], self.T_all[index], self.vlos_all[index]


class Forward(object):
    def __init__(self, gpu=0, checkpoint=None, readir=None, verbose=0):

        self.logger = logging.getLogger("neural")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []
        ch = logging.StreamHandler()        
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        # Is a GPU available?
        self.cuda = torch.cuda.is_available()
        self.gpu = gpu
        self.device = torch.device("cpu") #f"cuda:{self.gpu}" if self.cuda else "cpu")        

        if (checkpoint is None):
            if readir is None:
                raise ValueError('Not checkpoint or read directory selected')
            files = glob.glob(readir + '*.pth')
            self.checkpoint = sorted(files)[-1]
        else:
            self.checkpoint = checkpoint
            
        if (self.cuda):
            checkpoint = torch.load(self.checkpoint, map_location=lambda storage, loc: storage)
        else:
            checkpoint = torch.load(self.checkpoint, map_location=lambda storage, loc: storage)

        self.hyperparameters = checkpoint['hyperparameters']
        self.predict_model = graphnet.EncodeProcessDecode(**self.hyperparameters).to(self.device)
        self.predict_model.load_state_dict(checkpoint['state_dict'])

        self.predict_model.eval()

        if (verbose >= 1):
            npars = sum(p.numel() for p in self.predict_model.parameters() if p.requires_grad)
            tmp = self.checkpoint.split('/')
            self.logger.info(f'Using neural checkpoint {tmp[-1]} on {self.device} - N. parameters = {npars}')

    def predict(self, tau_all, ne_all, vturb_all, T_all, vlos_all):

        dset = Dataset(self.hyperparameters, tau_all, ne_all, vturb_all, T_all, vlos_all)

        self.loader = torch_geometric.loader.DataLoader(dset, batch_size=1, shuffle=False)
        
        self.pred_out = []

        with torch.no_grad():

            for data in self.loader:

                node = data.x
                edge_attr = data.edge_attr
                edge_index = data.edge_index
                u = data.u
                batch = data.batch

                node, edge_attr, edge_index = node.to(self.device), edge_attr.to(self.device), edge_index.to(self.device)
                u, batch = u.to(self.device), batch.to(self.device)

                out = self.predict_model(node, edge_attr, edge_index, u, batch)

                n = len(data.ptr) - 1
                for i in range(n):
                    left = data.ptr[i]
                    right = data.ptr[i+1]
                    self.pred_out.append(out[left:right, :].cpu().numpy()*5)

        return self.pred_out
