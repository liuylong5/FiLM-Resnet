import torch
import torch.nn as nn
import torch.nn.functional as F

class ResNetBlock(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(ResNetBlock, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, input_size)

    def forward(self, x, gamma=None, beta=None):
        residual = x
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        if gamma is not None and beta is not None:
            out = gamma * out + beta
        out += residual
        out = self.relu(out)
        return out


class FiLMLight(nn.Module):
    def __init__(self, cond_dim, feature_dim):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(cond_dim, 64),
            nn.ReLU(),
            nn.Linear(64, feature_dim * 2)
        )

    def forward(self, cond):
        gamma_beta = self.mlp(cond)
        gamma, beta = gamma_beta.chunk(2, dim=-1)
        return gamma, beta


class PreprocessLayer(nn.Module):
    def __init__(self, Xmu, Xstd, Dmu, Dstd, Pmin, Pmax):
        super(PreprocessLayer, self).__init__()
        self.register_buffer('Xmu', Xmu)
        self.register_buffer('Xstd', Xstd)
        self.register_buffer('Dmu', Dmu)
        self.register_buffer('Dstd', Dstd)
        self.register_buffer('Pmin', Pmin)
        self.register_buffer('Pmax', Pmax)
        self.register_buffer('lambda_val', torch.Tensor([
        0.1000, 0.1000, 0.1000, 0.5000, 0.1000, 0.1000, 0.1000, 0.2000, 0.1000,
        0.1000, 0.0500, 0.1000, 0.3000, 0.2000, 0.1000, 0.1000, 0.2000, 0.0500,
        0.0100, 0.1000, 0.0100, 0.1000, 0.1000, 0.3000, 0.1000, 0.0500, 0.1000,
        0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.0500, 0.0800, 0.1000,
        0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000,
        0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1500, 0.0500, 0.1000,
        0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.1000, 0.8000, 0.1000, 0.1000,
        0.1000, 0.1000, 0.1000]))

        no_change_idx = [48, 50, 53, 56, 57, 58, 59, 62, 63, 65]
        all_indices = set(range(66))
        self.change_indices = sorted(list(all_indices - set(no_change_idx)))
        self.lambda_val_init = self.lambda_val[self.change_indices]

    def box_cox_transform(self, x, lambda_val):
        lambda_val = lambda_val.to(x.device)
        return (torch.pow(x, lambda_val) - 1) / lambda_val

    def inverse_box_cox(self, x_bct, lambda_val):
        lambda_val = lambda_val.to(x_bct.device)
        return torch.pow(x_bct * lambda_val + 1, 1 / lambda_val)


    def forward(self, init_conc, params, final_conc, all_init_conc):

        x_bct = self.box_cox_transform(init_conc, self.lambda_val_init)
        x_norm = (x_bct - self.Xmu) / self.Xstd
        p_norm = (params - self.Pmin) / (self.Pmax - self.Pmin)
        bct_init_all = self.box_cox_transform(all_init_conc, self.lambda_val)
        bct_final = self.box_cox_transform(final_conc, self.lambda_val)

        D = bct_final - bct_init_all
        d_norm = (D - self.Dmu) / self.Dstd
        return x_norm, p_norm, d_norm


    def inverse_transform(self, d_pred_norm, x_norm):
        batch_size = d_pred_norm.shape[0]
        device = d_pred_norm.device

        D = d_pred_norm * self.Dstd + self.Dmu

        x_bct = x_norm * self.Xstd + self.Xmu
        init_conc_p = self.inverse_box_cox(x_bct, self.lambda_val_init)

        full_init = torch.full((batch_size, 66), 1e-19, device=device)
        full_init[:, self.change_indices] = init_conc_p
        bct_init_full = self.box_cox_transform(full_init, self.lambda_val)
        bct_final = bct_init_full + D
        final_conc = self.inverse_box_cox(bct_final, self.lambda_val)

        final_conc = torch.clamp(final_conc, min=1e-19) 
        delta_conc = final_conc - full_init

        return delta_conc, final_conc



class DNN(nn.Module):
    def __init__(self, input_size, param_size, output_size, norm_args, use_film=True):
        super(DNN, self).__init__()
        self.loss_fun = nn.MSELoss(reduction='none')
        # self.loss_fun = nn.HuberLoss(reduction='none', delta=0.1)
        self.preprocess = PreprocessLayer(*norm_args)

        self.use_film = use_film
        self.input_layer = nn.Linear(input_size + param_size + 2, 1024)
        self.relu = nn.GELU()

        self.res_blocks = nn.ModuleList([
            ResNetBlock(1024, 512) for _ in range(6)
        ])

        if self.use_film:
            self.film = FiLMLight(cond_dim=1, feature_dim=1024)  

        self.output_layer = nn.Linear(1024, output_size)

    def forward(self, init_conc, params, final_conc, all_init_conc, if_light, min_step, return_delta=False):
        x_norm, p_norm, d_norm = self.preprocess(init_conc, params, final_conc, all_init_conc)
        x = torch.cat((x_norm, p_norm, if_light, min_step), dim=1)
        out = self.input_layer(x)
        out = self.relu(out)

        if self.use_film:
            gamma, beta = self.film(if_light)
        else:
            gamma = beta = None

        for block in self.res_blocks:
            out = block(out, gamma, beta)

        d_pred_norm = self.output_layer(out)

        if return_delta:
            return self.preprocess.inverse_transform(d_pred_norm, x_norm)
        else:
            return d_pred_norm, d_norm

    def loss(self, d_pred_norm, d_norm):
        return self.loss_fun(d_pred_norm, d_norm).mean()
