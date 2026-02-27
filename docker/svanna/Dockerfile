FROM maven:3.8.5-openjdk-17-slim AS build

RUN apt-get update && apt-get install -y git gzip wget bash && rm -rf /var/lib/apt/lists/*

RUN git clone --branch v1.0.5 --depth 1 https://github.com/TheJacksonLaboratory/SvAnna.git /app/svAnna

WORKDIR /app/svAnna

RUN mvn package

FROM eclipse-temurin:17-jre-slim

RUN apt-get update && apt-get install -y gzip bash wget p7zip-full && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /app

RUN wget -O /app/2304_hg38.svanna.zip https://storage.googleapis.com/svanna/2304_hg38.svanna.zip \
    && 7z x /app/2304_hg38.svanna.zip -o/app/svanna_db_2304_hg38 -y \
    && rm /app/2304_hg38.svanna.zip

COPY --from=build /app/svAnna/svanna-cli/target/svanna-cli-*.jar /app/svanna-cli.jar

WORKDIR /app

CMD ["java", "-jar", "svanna-cli.jar", "prioritize", "-d", "/app/svanna_db_2304_hg38"]