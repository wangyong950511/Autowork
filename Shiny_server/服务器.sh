# 更新脚本到服务器
sudo cp qpcr_module.R protein_module.R primer_module.R app.R /srv/shiny-server/autorun/

# 查看日志
tail -f $(ls -t /var/log/shiny-server/myapp-*.log | head -n 1)

# 重启shiny-server
sudo systemctl restart shiny-server
