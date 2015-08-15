import matplotlib.pyplot as plt
import numpy
import sys

if len(sys.argv) < 3:
    print "Usage: "
    print "    python %s <filename> <label> [<filename> <label>]" % sys.argv[0]
    exit()

filenames = sys.argv[1::2]
labels    = sys.argv[2::2]

for i in range(len(filenames)):
    x,rho,u,p,mu,bulk,kap = numpy.loadtxt(filenames[i], skiprows=2, unpack=True)

    plt.figure(1)
    plt.plot(x,rho,label=labels[i])

    plt.figure(2)
    plt.plot(x,u,label=labels[i])

    plt.figure(3)
    plt.plot(x,p,label=labels[i])

    plt.figure(4)
    plt.plot(x,mu,label=labels[i])

    plt.figure(5)
    plt.plot(x,bulk,label=labels[i])

    plt.figure(6)
    plt.plot(x,kap,label=labels[i])

plt.figure(1)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$\\rho$', fontsize=20)
plt.legend(loc='lower left')
plt.tight_layout()

plt.figure(2)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$u$', fontsize=20)
plt.legend(loc='lower left')
plt.tight_layout()

plt.figure(3)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$p$', fontsize=20)
plt.legend(loc='lower left')
plt.tight_layout()

plt.figure(4)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$\\mu$', fontsize=20)
plt.legend(loc='upper left')
plt.tight_layout()

plt.figure(5)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$\\beta$', fontsize=20)
plt.legend(loc='upper left')
plt.tight_layout()

plt.figure(6)
plt.xlim((x[0],x[-1]))
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$\\kappa$', fontsize=20)
plt.legend(loc='upper left')
plt.tight_layout()

plt.show()
