FROM nginx:stable-alpine

COPY subpath.template /etc/nginx/templates/default.conf.template

RUN mkdir -p /var/log/nginx /var/cache/nginx && \
    touch /var/log/nginx/api.log && \
    touch /var/log/nginx/error.log && \
    chmod 666 /var/log/nginx/api.log && \
    chmod 666 /var/log/nginx/error.log

CMD ["nginx", "-g", "daemon off;"]
