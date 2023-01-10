def plot_plaq(plaquette):
    rs = plaquette['rs']
    Rs = plaquette['Rs']
    xs = rs[:,0]
    ys = rs[:,1]
    print(xs)
    print(ys)
    for j, R in enumerate(plaquette['Rs']):
        if j == 0:
            color='black'
        else:
            color='red'
        for i in range(plaquette['L']):
            plt.scatter([xs[i]+R[0]], [ys[i]+R[1]], color=color)
            plt.text(xs[i]+R[0], ys[i]+R[1], str(i))