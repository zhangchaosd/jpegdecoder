#include<cstdio>
#include<cstdlib>
using namespace std;

#define COMPONENTS 3
#define HUFFMAN_TABLES 4
#define HUFFMAN_HASH_NBITS 9
#define HUFFMAN_HASH_SIZE  (1UL<<HUFFMAN_HASH_NBITS)
#define be16_to_cpu(x) (((x)[0]<<8)|(x)[1])
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef int jmp_buf[16];

enum std_markers {
	DQT = 0xDB, /* Define Quantization Table */
	SOF = 0xC0, /* Start of Frame (size information) */
	DHT = 0xC4, /* Huffman Table */
	SOI = 0xD8, /* Start of Image */
	SOS = 0xDA, /* Start of Scan */
	RST = 0xD0, /* Reset Marker d0 -> .. */
	RST7 = 0xD7, /* Reset Marker .. -> d7 */
	EOI = 0xD9, /* End of Image */
	DRI = 0xDD, /* Define Restart Interval */
	APP0 = 0xE0,
};

struct huffman_table
{
	/* Fast look up table, using HUFFMAN_HASH_NBITS bits we can have directly the symbol,
	* if the symbol is <0, then we need to look into the tree table */
	short int lookup[HUFFMAN_HASH_SIZE];
	/* code size: give the number of bits of a symbol is encoded */
	unsigned char code_size[HUFFMAN_HASH_SIZE];
	/* some place to store value that is not encoded in the lookup table
	* FIXME: Calculate if 256 value is enough to store all values
	*/
	uint16_t slowtable[16 - HUFFMAN_HASH_NBITS][256];
};
struct component
{
	unsigned int Hfactor;
	unsigned int Vfactor;
	float *Q_table;		/* Pointer to the quantisation table to use */
	struct huffman_table *AC_table;
	struct huffman_table *DC_table;
	short int previous_DC;	/* Previous DC coefficient */
	short int DCT[64];		/* DCT coef */
	unsigned int cid;
};
struct jdec_private
{
	/* Public variables */
	uint8_t *components[COMPONENTS];
	unsigned int width, height;	/* Size of the image */
	unsigned int flags;
	/* Private variables */
	const unsigned char *stream_begin, *stream_end;
	unsigned int stream_length;

	const unsigned char *stream;	/* Pointer to the current stream */
	unsigned int reservoir, nbits_in_reservoir;

	struct component component_infos[COMPONENTS];
	float Q_tables[COMPONENTS][64];		/* quantization tables */
	struct huffman_table HTDC[HUFFMAN_TABLES];	/* DC huffman tables   */
	struct huffman_table HTAC[HUFFMAN_TABLES];	/* AC huffman tables   */
	int default_huffman_table_initialized;
	int restart_interval;
	int restarts_to_go;				/* MCUs left in this restart interval */
	int last_rst_marker_seen;			/* Rst marker is incremented each time */

										/* Temp space used after the IDCT to store each components */
	uint8_t Y[64 * 4], Cr[64], Cb[64];

	jmp_buf jump_state;
	/* Internal Pointer use for colorspace conversion, do not modify it !!! */
	uint8_t *plane[COMPONENTS];

};

int getfilesize(FILE *f)
{
	fseek(f, 0, SEEK_END);
	long pos = ftell(f);
	fseek(f, 0, SEEK_SET);
	return pos;
}

static int parse_JFIF(struct jdec_private *priv, const unsigned char *stream)
{
	int chuck_len;
	int marker;
	int sos_marker_found = 0;
	int dht_marker_found = 0;
	const unsigned char *next_chunck;

	/* Parse marker */
	while (!sos_marker_found)
	{
		if (*stream++ != 0xff)
			goto bogus_jpeg_format;
		/* Skip any padding ff byte (this is normal) */
		while (*stream == 0xff)
			stream++;

		marker = *stream++;
		chuck_len = be16_to_cpu(stream);
		next_chunck = stream + chuck_len;
		switch (marker)
		{
		case SOF:
			if (parse_SOF(priv, stream) < 0)
				return -1;
			break;
		case DQT:
			if (parse_DQT(priv, stream) < 0)
				return -1;
			break;
		case SOS:
			if (parse_SOS(priv, stream) < 0)
				return -1;
			sos_marker_found = 1;
			break;
		case DHT:
			if (parse_DHT(priv, stream) < 0)
				return -1;
			dht_marker_found = 1;
			break;
		case DRI:
			if (parse_DRI(priv, stream) < 0)
				return -1;
			break;
		default:
			trace("> Unknown marker %2.2x\n", marker);
			break;
		}

		stream = next_chunck;
	}

	if (!dht_marker_found) {
		trace("No Huffman table loaded, using the default one\n");
		build_default_huffman_tables(priv);
	}

#ifdef SANITY_CHECK
	if ((priv->component_infos[cY].Hfactor < priv->component_infos[cCb].Hfactor)
		|| (priv->component_infos[cY].Hfactor < priv->component_infos[cCr].Hfactor))
		error("Horizontal sampling factor for Y should be greater than horitontal sampling factor for Cb or Cr\n");
	if ((priv->component_infos[cY].Vfactor < priv->component_infos[cCb].Vfactor)
		|| (priv->component_infos[cY].Vfactor < priv->component_infos[cCr].Vfactor))
		error("Vertical sampling factor for Y should be greater than vertical sampling factor for Cb or Cr\n");
	if ((priv->component_infos[cCb].Hfactor != 1)
		|| (priv->component_infos[cCr].Hfactor != 1)
		|| (priv->component_infos[cCb].Vfactor != 1)
		|| (priv->component_infos[cCr].Vfactor != 1))
		error("Sampling other than 1x1 for Cr and Cb is not supported");
#endif

	return 0;
bogus_jpeg_format:
	trace("Bogus jpeg format\n");
	return -1;
}

int tinyjpeg_parse_header(struct jdec_private *priv, const unsigned char *buf, unsigned int size)
{
	int ret;

	/* Identify the file */
	if ((buf[0] != 0xFF) || (buf[1] != SOI))
		printf("Not a JPG file ?\n");

	priv->stream_begin = buf + 2;
	priv->stream_length = size - 2;
	priv->stream_end = priv->stream_begin + priv->stream_length;

	ret = parse_JFIF(priv, priv->stream_begin);

	return ret;
}

int main() {
	char url[] = "test.jpg";
	FILE *fp = fopen(url, "rb");
	if (!fp)
	{
		printf("jpg open failed.\n");
		return -1;
	}
	int length_of_file = getfilesize(fp);//获得文件长度
	unsigned char *buf;
	buf = (unsigned char *)malloc(length_of_file + 4);
	if (!buf)
	{
		printf("alloc buf memory failed.\n");
		return -1;
	}
	fread(buf, length_of_file, 1, fp);
	fclose(fp);
	struct jdec_private *jdec;
	jdec = (struct jdec_private *)calloc(1, sizeof(struct jdec_private));
	//jdec->dlg=(CSpecialVIJPGDlg *)lparam;?????????????


	printf("%d", 12);
	return 0;
}