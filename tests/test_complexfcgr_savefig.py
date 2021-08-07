from pathlib import Path
import random
from complexcgr import complexFCGR

BASEPATH = Path().parent.resolve()
def test_savefig():
    fcgr = complexFCGR(k=8)
    seq = "".join(random.choice("ACG") for _ in range(300_000))
    fig = fcgr(seq)
    fcgr.save(fig, path=BASEPATH.joinpath("img/ACG-complexCGR.png"))
