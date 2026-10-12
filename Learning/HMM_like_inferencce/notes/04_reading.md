# 阅读指南与原始来源

下面是阅读导航，讲义与代码为独立教学推导，非书籍逐章改写。按主题定位 Murphy 章节，避免不同版本页码差异。

1. [Bishop, PRML（作者/出版社提供的 PDF）](https://www.microsoft.com/en-us/research/wp-content/uploads/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf)：第 13 章 sequential data，先 HMM/forward-backward/EM，再 linear dynamical systems。VI 背景读第 10 章。输出：自己推导 γ,ξ，不只读公式。
2. [Murphy, Probabilistic Machine Learning: Advanced Topics](https://probml.github.io/book2)：查 state-space models、variational inference、variational autoencoders；与 Bishop 对照符号。输出：区分 discrete HMM、linear Gaussian SSM、nonlinear/non-Gaussian SSM。
3. [Ghahramani & Jordan, Factorial Hidden Markov Models, 1997](https://mlg.eng.cam.ac.uk/pub/pdf/GhaJor97a.pdf)：重点独立 latent chains 与 coupled posterior，以及 approximate inference。输出：两个 emitter 的 graphical model。
4. [Kingma & Welling, Auto-Encoding Variational Bayes](https://arxiv.org/abs/1312.6114)：ELBO、recognition model、reparameterization。输出：Gaussian KL 推导与计算图。
5. [PyTorch 官方 autograd 教程](https://docs.pytorch.org/tutorials/beginner/introyt/autogradyt_tutorial.html)：requires_grad、计算图、grad、detach。输出：固定样本 score 与 reparameterized gradient 的区别。
6. [Johnson et al., 2016, SVAE](https://papers.nips.cc/paper_files/paper/2016/hash/7d6044e95a16761171b130dcb476a43e-Abstract.html)：先读 graphical model + recognition potentials 的结构，再看 objective。输出：网络负责哪一步、图模型负责哪一步。
7. [Krishnan et al., Deep Kalman Filters, 2015](https://arxiv.org/abs/1511.05121) 和 [Structured Inference Networks for Nonlinear State Space Models, 2017](https://arxiv.org/abs/1609.09869)：区分生成结构与 posterior approximation，关注 future observations 的使用。
8. [Cranmer, Brehmer & Louppe, PNAS 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7720103/)：辨别 posterior/likelihood/ratio estimation。输出：何时已有可求 likelihood，何时确实需要 SBI。

CRLB 学习以 `03_information.md` 和 notebooks 05–06 的可核对推导为主；进一步用统计教材的 likelihood、Fisher information、regularity、nuisance parameters 章节补理论，不从神经网络训练 loss 反推 Fisher。
