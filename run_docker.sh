running_image=$(sudo docker images | grep 'dis_prs_calc' | wc -l)
if [[ $running_image -lt 1 ]]; then
    sudo docker build -t dis_prs_calc .
fi
running_docker=$(sudo docker container ls | grep 'dis_prs_calc_cont' | wc -l)
stopped_docker=$(sudo docker container ls -a | grep 'dis_prs_calc_cont' | wc -l)

if [[ $running_docker -lt 1 && $stopped_docker -lt 1 ]]; then
    sudo docker run --name dis_prs_calc_cont -t -d dis_prs_calc 
elif [[ $running_docker -lt 1 && $stopped_docker -ge 1 ]]; then
    sudo docker start dis_prs_calc_cont
fi


sudo docker exec dis_prs_calc_cont python3 manage.py runserver 0.0.0.0:8000 & echo "Your container is running. You can log in as user (name: testuser, pw: 1234) or admin (name: admin, pw: 1234) at: " $(sudo docker inspect --format '{{ .NetworkSettings.IPAddress }}' dis_prs_calc_cont):8000/dis_calc/login.html
