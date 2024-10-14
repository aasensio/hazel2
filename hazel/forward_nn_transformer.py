import numpy as np
import glob
import torch
import torch.nn as nn
from sklearn import neighbors
import logging


class PositionalEncoding(nn.Module):
    def __init__(self, d_emb, norm=10000.0):
        """
        Inputs
            d_model - Hidden dimensionality.
        """
        super().__init__()
        self.d_emb = d_emb
        self.norm = norm

    def forward(self, t):
        pe = torch.zeros(t.shape[0], t.shape[1], self.d_emb).to(t.device)  # (B, T, D)
        div_term = torch.exp(torch.arange(0, self.d_emb, 2).float() * (-np.log(self.norm) / self.d_emb))[None, None, :].to(t.device)  # (1, 1, D / 2)
        t = t.unsqueeze(2)  # (B, 1, T)
        pe[:, :, 0::2] = torch.sin(t * div_term)  # (B, T, D / 2)
        pe[:, :, 1::2] = torch.cos(t * div_term)  # (B, T, D / 2)
        return pe  # (B, T, D)
    
class TransformerModel(nn.Module):

    def __init__(self, ninp, nemb, nout, nhead, nhid, nlayers, dropout=0.1, norm=1000.0):
        """
        Transformer model for sequence to sequence learning

        Args:
            ninp (_type_): input size
            nemb (_type_): embedding size
            nout (_type_): output size
            nhead (_type_): number of heads
            nhid (_type_): hidden layer size in feed forward network
            nlayers (_type_): number of layers
            dropout (float, optional): dropout probability. Defaults to 0.5.
        """
        super(TransformerModel, self).__init__()

        self.model_type = 'Transformer'
        
        self.encoder = nn.Linear(ninp, nemb)

        self.pos_encoder = PositionalEncoding(nemb, norm)
        
        encoder_layers = nn.TransformerEncoderLayer(nemb, nhead, nhid, dropout, norm_first=True, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, nlayers, enable_nested_tensor=False)
        
        self.nemb = nemb
        self.decoder = nn.Linear(nemb, nout)

        self.init_weights()

    def init_weights(self):

        # Since TraansformerEncoder inputs a TransformerEncoderLayer, all layers will use exactly the same initialization
        # We undo this here
        for name, param in self.named_parameters():
            if 'weight' in name and param.data.dim() == 2:
                nn.init.kaiming_uniform_(param)

    def forward(self, src, tau, src_mask):
                        
        # Get tau embedding
        tau_emb = self.pos_encoder(tau)

        # Embed the input sequence into the embedding space and add the tau embedding
        x = self.encoder(src) + tau_emb
            
        # Apply the transformer encoder
        x = self.transformer_encoder(x, src_key_padding_mask=src_mask)

        # Apply the decoder to the output space
        x = self.decoder(x)

        output = (~src_mask).float()[:, :, None] * x

        return output


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
            
        checkpoint = torch.load(self.checkpoint, map_location=lambda storage, loc: storage, weights_only=False)
        
        self.hyperparameters = checkpoint['hyperparameters']
        self.predict_model = TransformerModel(ninp=self.hyperparameters['transformer']['n_input'],
                                            nemb=self.hyperparameters['transformer']['n_embedding'], 
                                            nout=self.hyperparameters['transformer']['n_output'], 
                                            nhead=self.hyperparameters['transformer']['n_heads'], 
                                            nhid=self.hyperparameters['transformer']['n_hidden'], 
                                            nlayers=self.hyperparameters['transformer']['n_layers'],
                                            norm=self.hyperparameters['transformer']['norm'],
                                            dropout=self.hyperparameters['transformer']['dropout']).to(self.device)        
        self.predict_model.load_state_dict(checkpoint['state_dict'])

        self.predict_model.eval()

        if (verbose >= 1):
            npars = sum(p.numel() for p in self.predict_model.parameters() if p.requires_grad)
            tmp = self.checkpoint.split('/')
            self.logger.info(f'    * Using neural checkpoint {tmp[-1]} on {self.device} - N. parameters = {npars}')

    def predict(self, tau_all, ne_all, vturb_all, T_all, vlos_all):        

        tau = (np.log10(tau_all.astype('float32')) + 10.0) * 10.0
        vturb = vturb_all.astype('float32') / 1e3 - 6.0
        vlos = vlos_all.astype('float32') / 1e3
        T = np.log10(T_all.astype('float32')) - 3.8
        ne = np.log10(ne_all.astype('float32')) - 16.0

        pars = np.concatenate([vturb[None, :], vlos[None, :], T[None, :], ne[None, :]], axis=0).T
        mask = np.zeros(len(tau)).astype('bool')

        pars = torch.tensor(pars, dtype=torch.float32).to(self.device)
        tau = torch.tensor(tau, dtype=torch.float32).to(self.device)
        mask = torch.tensor(mask, dtype=torch.bool).to(self.device)
        
        with torch.no_grad():
            self.pred_out = self.predict_model(pars[None, ...], tau[None, ...], mask[None, ...])
        
        return self.pred_out[0, ...].cpu().numpy()
