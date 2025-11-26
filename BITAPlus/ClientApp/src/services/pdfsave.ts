import { Injectable } from '@angular/core';
import html2canvas from 'html2canvas';
import jsPDF from 'jspdf';


@Injectable({
  providedIn: 'root'
})
export class PdfSave {
  async exportToPdf(element: HTMLElement | null, filename: string = 'export.pdf', marginMM: number = 10, targetDPI: number = 300): Promise<void> {
    try {
      if (!element) throw new Error('No element provided');

      // Capture element as canvas
      const canvas = await html2canvas(element, {
        scale: targetDPI / 96,
        useCORS: true,
        allowTaint: true,
        backgroundColor: '#ffffff',
      });

      // Get canvas dimensions
      const imgWidth = canvas.width;
      const imgHeight = canvas.height;
      const aspect = imgHeight / imgWidth;

      // Create PDF
      const doc = new jsPDF();
      const pageW = doc.internal.pageSize.getWidth();
      const pageH = doc.internal.pageSize.getHeight();
      const desiredW_MM = pageW - marginMM * 2;
      const desiredH_MM = desiredW_MM * aspect;

      // Convert canvas to image and add to PDF
      const imgData = canvas.toDataURL('image/png');
      doc.addImage(imgData, 'PNG', marginMM, marginMM, desiredW_MM, desiredH_MM);

      // Handle multi-page content
      let heightLeft = desiredH_MM - (pageH - marginMM * 2);
      let position = pageH - marginMM;
      while (heightLeft > 0) {
        doc.addPage();
        doc.addImage(imgData, 'PNG', marginMM, position, desiredW_MM, desiredH_MM);
        heightLeft -= (pageH - marginMM * 2);
        position = pageH - marginMM;
      }

      doc.save(filename);
      console.log('PDF generated successfully');
    } catch (error) {
      console.error('Error generating PDF:', error);
      throw error;
    }
  }
}