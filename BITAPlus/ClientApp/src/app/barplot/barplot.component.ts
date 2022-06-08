import { Component, Input, ElementRef, OnChanges, OnInit, Output, EventEmitter } from '@angular/core';
import * as d3 from 'd3';
export const convHackForMinMax=(a:number|undefined)=>a!=undefined?a:0;

@Component({
  selector: 'app-barplot',
  templateUrl: './barplot.component.html',
  styleUrls: ['./barplot.component.css']
})
export class BarplotComponent implements OnChanges, OnInit {
  @Input() DATA: Array<number> = [];
  @Input() editdata: Array<number> = [];
  @Output() editChange: EventEmitter<Array<number>> = new EventEmitter<Array<number>>();
  dataToChange: Array<number> = [];
  @Input() width = 500;
  @Input() height = 500;
  @Input() title = '';
  @Input() editvalues = false;
  abshack = Math.abs;
  scaleX = d3.scaleLinear();
  scaleY = d3.scaleLinear();
  rimX = 0.15;
  rimY = 0.1;
  format = (n: number) => d3.format('2.4f')(n);
  formatL = (n: number) => d3.format('2.1f')(n);
  translatehack = (x: number, y: number, r = 0) => `translate(${x},${y}) rotate(${r})`;
  constructor(private element: ElementRef) { }

  ngOnChanges() {
    this.dataToChange = this.DATA.map(d => d);
    console.log(this.DATA);
    setTimeout(() => {
      this.update();
    }, 100);
  }
  ngOnInit() {
    this.dataToChange = this.DATA.map(d => d);
  }
  info(e: MouseEvent, x: number, y: number, inout = false) {
    const tip = d3.select(this.element.nativeElement).select('div.mainTip');
    console.log(tip);
    const Torg = ((d3.select('body').node() as HTMLElement)
    /*.parentNode.parentNode.parentNode.parentNode.parentNode.parentNode as HTMLElement*/)
      .getBoundingClientRect().left;
    const origin = (d3
      .select('app-optimise')
      .node() as HTMLElement).getBoundingClientRect(); // Try to get position correct when the picture has scrollbars.
    const here = d3.select(e.target as HTMLElement & EventTarget);
    if (inout) {
      here.style('opacity', 0.5);
      tip // The tooltip
        .style('left', `${e.clientX - 60 - 0 * origin.left + 0 * this.width / this.DATA.length}px`)
        .style('top', `${e.clientY - origin.top}px`)
        .style('opacity', 1)
        .style('display', 'inline-block')
        .html(`x:${x} ${this.format(y)}`);
      const wwww = (tip.node() as HTMLElement).getBoundingClientRect();
      console.log(wwww.width, origin.left - Torg);
      tip.style('left', `${e.clientX - wwww.width * 0.5 - Torg}px`);
    } else {
      here.style('opacity', 1);
      tip.style('opacity', 0).style('display', 'none');
    }
  }
  getnewChange = (index: number) => {
    const newtext = +(d3.select(this.element.nativeElement).select('input.newdat').node() as HTMLInputElement).value;
    const textv = d3.select(this.element.nativeElement).selectAll('text.editbox').nodes()[index] as SVGTextElement;
    this.dataToChange[index] = newtext;
    console.log(this.dataToChange);
    this.editdata = this.dataToChange;
    this.update();
  }

  update() {
    this.scaleX
      .domain([0, this.dataToChange.length])
      .range([this.width * this.rimX, this.width * (1 - 0.5 * this.rimX)]);
    const useforscale = this.editvalues ? this.dataToChange : this.DATA;
    this.scaleY
      .domain([Math.min(convHackForMinMax(d3.min(useforscale)), 0), convHackForMinMax(d3.max(useforscale))])
      .range([this.height * (1 - this.rimY), this.height * this.rimY]);
    d3.select(this.element.nativeElement).selectAll('rect.ongraph').data(this.DATA)
      .transition()
      .duration(2000)
      .attrTween('y', d => t => {
        return `${t * this.scaleY(d >= 0 ? d : d * (1 - t))}`;
      })
      .attrTween('height', (d, i) => t => {
        return `${this.abshack(this.scaleY(d * t) - this.scaleY(d * (1 - t)))}`;
      });
    this.editChange.emit(this.editdata);
  }
}