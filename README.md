# 2024年4月14日第一次建立仓库
## 解决ssh代理问题
Host https://github.com    
	ProxyCommand connect -S 127.0.0.1:7890 %h %p
